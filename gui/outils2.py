from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure
from PySide6.QtGui import QFontDatabase
from PySide6.QtWidgets import (QCheckBox, QComboBox, QFileDialog, QHBoxLayout, QLabel,
                               QMessageBox, QPlainTextEdit, QPushButton, QSpinBox,
                               QVBoxLayout, QWidget)

from core.alignement import type_sequence
from core.annotation import chercher_orfs, profil_hydropathie, regions_hydrophobes, traduire
from core.fasta import lire_fasta
from core.msa import aligner_multiple, conservation, trouver_programme, vers_fasta
from gui.travailleur import Travailleur


def _zone_texte():
    zone = QPlainTextEdit()
    zone.setReadOnly(True)
    zone.setFont(QFontDatabase.systemFont(QFontDatabase.FixedFont))
    return zone


# ---------------------------------------------------------------- traduction / ORF
def page_orf(dialogue):
    page = QWidget()
    layout = QVBoxLayout(page)
    ligne = QHBoxLayout()
    cible = QComboBox()
    cible.addItems(["A", "B"])
    longueur = QSpinBox()
    longueur.setRange(1, 5000)
    longueur.setValue(30)
    deux = QCheckBox("6 cadres")
    deux.setChecked(True)
    cadre = QComboBox()
    cadre.addItems(["1", "2", "3", "-1", "-2", "-3"])
    btn_orf = QPushButton("Chercher les ORF")
    btn_tr = QPushButton("Traduire le cadre")
    for w in (QLabel("Séquence :"), cible, QLabel("ORF ≥"), longueur, QLabel("aa"), deux,
              btn_orf, cadre, btn_tr):
        ligne.addWidget(w)
    layout.addLayout(ligne)
    sortie = _zone_texte()
    layout.addWidget(sortie)

    def chercher():
        entree = dialogue._sequence(cible.currentText())
        if entree is None:
            return
        ident, seq = entree
        if type_sequence(seq) != "nucleotide":
            QMessageBox.information(dialogue, "ORF", "Cette séquence n'est pas de l'ADN.")
            return
        orfs = chercher_orfs(seq, longueur.value(), deux_brins=deux.isChecked())
        lignes = [f"{ident} : {len(orfs)} ORF complet(s) (ATG → stop)", ""]
        for o in orfs:
            prot = o["proteine"] if len(o["proteine"]) <= 60 else o["proteine"][:57] + "..."
            lignes.append(f"{o['brin']} cadre {o['cadre']:>2}  {o['debut']}–{o['fin']}  "
                          f"{o['longueur_aa']} aa  {prot}")
        sortie.setPlainText("\n".join(lignes))

    def traduire_cadre():
        entree = dialogue._sequence(cible.currentText())
        if entree is None:
            return
        ident, seq = entree
        if type_sequence(seq) != "nucleotide":
            QMessageBox.information(dialogue, "Traduction", "Cette séquence n'est pas de l'ADN.")
            return
        prot = traduire(seq, int(cadre.currentText()))
        blocs = [prot[i:i + 60] for i in range(0, len(prot), 60)]
        sortie.setPlainText(f"{ident}, cadre {cadre.currentText()} ({len(prot)} aa)\n\n"
                            + "\n".join(blocs))

    btn_orf.clicked.connect(chercher)
    btn_tr.clicked.connect(traduire_cadre)
    return page


# ---------------------------------------------------------------- hydropathie
def page_hydro(dialogue):
    page = QWidget()
    layout = QVBoxLayout(page)
    ligne = QHBoxLayout()
    cible = QComboBox()
    cible.addItems(["A", "B"])
    fenetre = QSpinBox()
    fenetre.setRange(1, 99)
    fenetre.setSingleStep(2)
    fenetre.setValue(19)
    bouton = QPushButton("Tracer (Kyte-Doolittle)")
    for w in (QLabel("Protéine :"), cible, QLabel("Fenêtre :"), fenetre, bouton):
        ligne.addWidget(w)
    ligne.addStretch(1)
    layout.addLayout(ligne)
    figure = Figure(figsize=(5, 3))
    canvas = FigureCanvasQTAgg(figure)
    layout.addWidget(canvas, 1)
    sortie = _zone_texte()
    sortie.setMaximumHeight(110)
    layout.addWidget(sortie)

    def tracer():
        entree = dialogue._sequence(cible.currentText())
        if entree is None:
            return
        ident, seq = entree
        if type_sequence(seq) == "nucleotide":
            QMessageBox.information(dialogue, "Hydropathie",
                                    "Séquence nucléique : traduis-la d'abord (onglet ORF).")
            return
        try:
            profil = profil_hydropathie(seq, fenetre.value())
            regions = regions_hydrophobes(seq, fenetre.value())
        except ValueError as e:
            QMessageBox.warning(dialogue, "Hydropathie", str(e))
            return
        figure.clear()
        ax = figure.add_subplot(111)
        ax.plot([p for p, _ in profil], [s for _, s in profil], linewidth=1.2)
        ax.axhline(1.6, linestyle="--", linewidth=0.8, color="#d6453d")
        for debut, fin in regions:
            ax.axvspan(debut, fin, alpha=0.2, color="#d98a00")
        ax.set_xlabel("position")
        ax.set_ylabel("Kyte-Doolittle")
        ax.set_title(ident)
        figure.tight_layout()
        canvas.draw()
        txt = [f"Fenêtre {fenetre.value()} ; seuil 1,6 (hélice transmembranaire : fenêtre 19)"]
        txt += [f"région hydrophobe {d}–{f} ({f - d + 1} aa)" for d, f in regions]
        if not regions:
            txt.append("aucune région au-dessus du seuil")
        sortie.setPlainText("\n".join(txt))

    bouton.clicked.connect(tracer)
    return page


# ---------------------------------------------------------------- alignement multiple
class PageMSA(QWidget):
    def __init__(self, dialogue):
        super().__init__()
        self.dialogue = dialogue
        self.entrees = []
        self.resultat = None
        self.tache = None
        layout = QVBoxLayout(self)
        ligne = QHBoxLayout()
        btn_charger = QPushButton("Charger un multi-FASTA…")
        btn_charger.clicked.connect(self.charger)
        self.avec_ab = QCheckBox("Ajouter A et B")
        self.btn_aligner = QPushButton("Aligner")
        self.btn_aligner.clicked.connect(self.aligner)
        self.btn_export = QPushButton("Exporter FASTA…")
        self.btn_export.clicked.connect(self.exporter)
        for w in (btn_charger, self.avec_ab, self.btn_aligner, self.btn_export):
            ligne.addWidget(w)
        ligne.addStretch(1)
        layout.addLayout(ligne)
        programme = trouver_programme()
        self.etat = QLabel("Fichier : (rien chargé) | programme : "
                           + (programme if programme else
                              "introuvable — sudo apt install clustalo"))
        layout.addWidget(self.etat)
        self.sortie = _zone_texte()
        layout.addWidget(self.sortie)

    def charger(self):
        chemin, _ = QFileDialog.getOpenFileName(self, "Multi-FASTA", "", "Tous (*)")
        if not chemin:
            return
        try:
            self.entrees = lire_fasta(chemin)
        except Exception as e:
            QMessageBox.warning(self, "Erreur", f"Lecture impossible : {e}")
            return
        self.etat.setText(f"Fichier : {len(self.entrees)} séquence(s) | "
                          f"programme : {trouver_programme() or 'introuvable'}")

    def _toutes(self):
        seqs = list(self.entrees)
        if self.avec_ab.isChecked():
            for cle in ("A", "B"):
                e = self.dialogue.fenetre.seqs.get(cle)
                if e:
                    seqs.append(e)
        return seqs

    def aligner(self):
        if self.tache is not None and self.tache.isRunning():
            return
        seqs = self._toutes()
        if len(seqs) < 2:
            QMessageBox.information(self, "Alignement multiple",
                                    "Charge un multi-FASTA (ou coche « Ajouter A et B »).")
            return
        self.btn_aligner.setEnabled(False)
        self.sortie.setPlainText("Calcul en cours…")
        self.tache = Travailleur(aligner_multiple, seqs)
        self.tache.termine.connect(self._fin)
        self.tache.echec.connect(self._echec)
        self.tache.start()

    def _echec(self, message):
        self.tache.wait()
        self.btn_aligner.setEnabled(True)
        self.sortie.setPlainText(f"Erreur : {message}")

    def _fin(self, alignement):
        self.tache.wait()
        self.btn_aligner.setEnabled(True)
        self.resultat = alignement
        cons = conservation(alignement)
        largeur_nom = min(max(len(i) for i, _ in alignement), 20)
        n = len(cons)
        lignes = [f"{len(alignement)} séquences, {n} colonnes, "
                  f"{cons.count('*')} colonnes conservées", ""]
        for debut in range(0, n, 60):
            for ident, seq in alignement:
                lignes.append(f"{ident[:largeur_nom]:<{largeur_nom}}  {seq[debut:debut + 60]}")
            lignes.append(f"{'':<{largeur_nom}}  {cons[debut:debut + 60]}")
            lignes.append("")
        self.sortie.setPlainText("\n".join(lignes))

    def exporter(self):
        if not self.resultat:
            QMessageBox.information(self, "Export", "Aucun alignement à exporter.")
            return
        chemin, _ = QFileDialog.getSaveFileName(self, "Exporter l'alignement",
                                                "", "FASTA (*.fasta *.fa)")
        if not chemin:
            return
        if "." not in chemin.rsplit("/", 1)[-1]:
            chemin += ".fasta"
        with open(chemin, "w", encoding="utf-8") as f:
            f.write(vers_fasta(self.resultat))
        self.dialogue.fenetre.statusBar().showMessage(f"Alignement exporté : {chemin}", 5000)

    def attendre(self):
        if self.tache is not None and self.tache.isRunning():
            self.tache.wait(10000)


def page_msa(dialogue):
    page = PageMSA(dialogue)
    dialogue._page_msa = page
    return page
