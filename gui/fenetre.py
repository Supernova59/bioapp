from PySide6.QtWidgets import (QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
                               QPushButton, QLabel, QFileDialog, QMessageBox, QSpinBox,
                               QComboBox, QPlainTextEdit, QLineEdit, QDoubleSpinBox)
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure

from core.fasta import lire_fasta
from core.dotplot import dotplot, tracer_dotplot
from core.alignement import aligner_complet, type_sequence
from core.ncbi import telecharger_fasta
from gui.exports import ExportMixin
from gui.outils import OutilsMixin
from core.sequences import complement_inverse
from PySide6.QtWidgets import QCheckBox
from PySide6.QtCore import QSettings
from gui.travailleur import Travailleur
from PySide6.QtGui import QFontDatabase

MAX_CASES = 25_000_000  # garde-fou mémoire (~75 Mo)


class Fenetre(OutilsMixin, ExportMixin, QMainWindow):
    def __init__(self, presence=None):
        super().__init__()
        self.presence = presence
        self.tache = None
        self.setWindowTitle("bioapp")
        self.seqs = {"A": None, "B": None}  # (identifiant, séquence)
        self.labels = {}

        centre = QWidget()
        self.setCentralWidget(centre)
        layout = QVBoxLayout(centre)

        for cle in ("A", "B"):
            ligne = QHBoxLayout()
            bouton = QPushButton(f"Charger FASTA {cle}")
            bouton.clicked.connect(lambda _=False, c=cle: self.charger(c))
            self.labels[cle] = QLabel("(rien chargé)")
            ligne.addWidget(bouton)
            ligne.addWidget(self.labels[cle], 1)
            layout.addLayout(ligne)

        ligne_ncbi = QHBoxLayout()
        self.accession = QLineEdit()
        self.accession.setPlaceholderText("accession NCBI (ex. NM_000546)")
        self.email = QLineEdit(QSettings("bioapp", "bioapp").value("email", ""))
        self.email.setPlaceholderText("e-mail (requis par NCBI)")
        self.cible = QComboBox()
        self.cible.addItems(["A", "B"])
        self.btn_ncbi = QPushButton("Télécharger")
        self.btn_ncbi.clicked.connect(self.telecharger)
        for w in (self.accession, self.email, self.cible, self.btn_ncbi):
            ligne_ncbi.addWidget(w)
        layout.addLayout(ligne_ncbi)

        options = QHBoxLayout()
        options.addWidget(QLabel("Fenêtre :"))
        self.spin = QSpinBox()
        self.spin.setRange(1, 200)
        self.spin.setValue(8)
        options.addWidget(self.spin)
        self.btn_dot = QPushButton("Dotplot")
        self.btn_dot.setEnabled(False)
        self.btn_dot.clicked.connect(self.tracer)
        options.addWidget(self.btn_dot)
        self.chk_rc = QCheckBox("B en complément inverse")
        options.addWidget(self.chk_rc)
        self.mode = QComboBox()
        self.mode.addItems(["global", "semi-global", "local"])
        options.addWidget(self.mode)
        self.btn_align = QPushButton("Aligner")
        self.btn_align.setEnabled(False)
        self.btn_align.clicked.connect(self.lancer_alignement)
        options.addWidget(self.btn_align)
        options.addStretch(1)
        layout.addLayout(options)

        bareme = QHBoxLayout()
        self.params = {}
        for cle, nom, defaut, mini, maxi in (
                ("match", "Match", 1, 0, 20),
                ("mismatch", "Mismatch", -1, -20, 0),
                ("ouverture", "Gap (ouverture)", -2, -20, 0),
                ("extension", "Gap (extension)", -1, -20, 0)):
            bareme.addWidget(QLabel(nom + " :"))
            spin = QDoubleSpinBox()
            spin.setDecimals(1)
            spin.setSingleStep(0.5)
            spin.setRange(mini, maxi)
            spin.setValue(defaut)
            self.params[cle] = spin
            bareme.addWidget(spin)
        bareme.addStretch(1)
        layout.addLayout(bareme)

        export = QHBoxLayout()
        export.addWidget(QLabel("Exporter :"))
        for texte, action in (("Dotplot (PNG/SVG)", self.exporter_dotplot),
                              ("Alignement (texte)", self.exporter_texte),
                              ("Alignement (FASTA)", self.exporter_fasta)):
            bouton = QPushButton(texte)
            bouton.clicked.connect(action)
            export.addWidget(bouton)
        btn_outils = QPushButton("Outils…")
        btn_outils.clicked.connect(self.ouvrir_outils)
        export.addWidget(btn_outils)
        export.addStretch(1)
        layout.addLayout(export)

        self.figure = Figure(figsize=(5, 5))
        self.canvas = FigureCanvasQTAgg(self.figure)
        layout.addWidget(self.canvas, 1)

        self.sortie = QPlainTextEdit()
        self.sortie.setReadOnly(True)
        self.sortie.setFont(QFontDatabase.systemFont(QFontDatabase.FixedFont))
        self.sortie.setMaximumHeight(180)
        layout.addWidget(self.sortie)

    def charger(self, cle):
        chemin, _ = QFileDialog.getOpenFileName(
            self, "Ouvrir un FASTA", "", "FASTA (*.fasta *.fa *.fna *.faa);;Tous (*)")
        if not chemin:
            return
        self.charger_chemin(cle, chemin)

    def charger_chemin(self, cle, chemin):
        try:
            entrees = lire_fasta(chemin)
            if not entrees:
                raise ValueError("aucune séquence trouvée")
        except Exception as e:
            QMessageBox.warning(self, "Erreur", f"Lecture impossible : {e}")
            return
        self.seqs[cle] = entrees[0]  # première séquence du fichier
        ident, seq = entrees[0]
        suite = f" (+{len(entrees) - 1} ignorées)" if len(entrees) > 1 else ""
        unite = "pb" if type_sequence(seq) == "nucleotide" else "aa"
        self.labels[cle].setText(f"{ident} — {len(seq)} {unite}{suite}")
        pret = all(self.seqs.values())
        self.btn_dot.setEnabled(pret)
        self.btn_align.setEnabled(pret)

    def statut(self, texte):
        self.statusBar().showMessage(texte)
        if self.presence is not None:
            self.presence.maj(texte)

    def occupe(self, oui):
        for w in (self.btn_dot, self.btn_align, self.btn_ncbi):
            w.setEnabled(not oui)
        if not oui:
            pret = all(self.seqs.values())
            self.btn_dot.setEnabled(pret)
            self.btn_align.setEnabled(pret)

    def lancer(self, message, fonction, apres, *args):
        if self.tache is not None and self.tache.isRunning():
            return
        self.occupe(True)
        self.statut(message)
        self._apres = apres
        self.tache = Travailleur(fonction, *args)
        self.tache.termine.connect(self._fin)
        self.tache.echec.connect(self._echec)
        self.tache.start()

    def _fin(self, resultat):
        self.tache.wait()
        self.occupe(False)
        self.statut("Analyse de séquences")
        self._apres(resultat)

    def _echec(self, message):
        self.tache.wait()
        self.occupe(False)
        self.statut("Analyse de séquences")
        QMessageBox.warning(self, "Erreur", message)

    def closeEvent(self, event):
        if self.tache is not None and self.tache.isRunning():
            self.tache.wait(10000)
        event.accept()

    def tracer(self):
        (id_a, a), (id_b, b) = self.seqs["A"], self.seqs["B"]
        if len(a) * len(b) > MAX_CASES:
            QMessageBox.warning(self, "Trop grand",
                                "Séquences trop longues pour un dotplot complet.")
            return
        if self.chk_rc.isChecked():
            if type_sequence(b) != "nucleotide":
                QMessageBox.information(self, "Dotplot",
                                        "Le complément inverse ne s'applique qu'à l'ADN.")
                return
            b = complement_inverse(b)
            id_b = f"{id_b} (compl. inverse)"
        self.noms_dot = (id_a, id_b)
        self.lancer(f"Dotplot de {id_a} vs {id_b}", dotplot,
                    self._afficher_dotplot, a, b, self.spin.value())

    def _afficher_dotplot(self, matrice):
        self.figure.clear()
        tracer_dotplot(matrice, self.figure.add_subplot(111), *self.noms_dot)
        self.canvas.draw()

    def lancer_alignement(self):
        (id_a, a), (id_b, b) = self.seqs["A"], self.seqs["B"]
        if len(a) * len(b) > MAX_CASES:
            QMessageBox.warning(self, "Trop grand",
                                "Séquences trop longues pour cet alignement.")
            return
        mode = self.mode.currentText()
        self.mode_calcule = mode
        self.ids_align = (id_a, id_b)
        p = {k: v.value() for k, v in self.params.items()}
        if "proteine" in (type_sequence(a), type_sequence(b)):
            args = (mode, 1, -1, -10, -1, "BLOSUM62")
            self.info_bareme = "protéines : BLOSUM62, gap -10/-1"
        else:
            args = (mode, p["match"], p["mismatch"], p["ouverture"], p["extension"], None)
            self.info_bareme = (f"ADN : match {p['match']}, mismatch {p['mismatch']}, "
                                f"gap {p['ouverture']}/{p['extension']}")
        self.lancer(f"Alignement {mode} de {id_a} vs {id_b}", aligner_complet,
                    self._afficher_alignement, a.upper(), b.upper(), *args)

    def _afficher_alignement(self, resultat):
        score, la, lb, (da, fa), (db, fb) = resultat
        self.dernier = (self.ids_align, la, lb)
        milieu = "".join("|" if x == y else " " for x, y in zip(la, lb))
        ident = milieu.count("|")
        texte = [f"Mode : {self.mode_calcule} | score : {score} | identité : {ident}/{len(la)}",
                 self.info_bareme]
        if self.mode_calcule == "local" and la:
            texte.append(f"Région alignée : A {da + 1}–{fa} | B {db + 1}–{fb}")
        texte.append("")
        for i in range(0, len(la), 60):
            texte += [la[i:i + 60], milieu[i:i + 60], lb[i:i + 60], ""]
        self.sortie.setPlainText("\n".join(texte))

    def telecharger(self):
        email = self.email.text().strip()
        QSettings("bioapp", "bioapp").setValue("email", email)
        accession = self.accession.text().strip()
        self.cible_tel = self.cible.currentText()
        self.lancer(f"Téléchargement de {accession}", telecharger_fasta,
                    self._fichier_telecharge, accession, email)

    def _fichier_telecharge(self, chemin):
        self.charger_chemin(self.cible_tel, str(chemin))
