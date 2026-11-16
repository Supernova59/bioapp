from PySide6.QtWidgets import (QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
                               QPushButton, QLabel, QFileDialog, QMessageBox, QSpinBox,
                               QComboBox, QPlainTextEdit, QLineEdit)
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure

from core.fasta import lire_fasta
from core.dotplot import dotplot, tracer_dotplot
from core.alignement import aligner
from core.ncbi import telecharger_fasta
from PySide6.QtCore import QSettings
from PySide6.QtGui import QFontDatabase

MAX_CASES = 25_000_000  # garde-fou mémoire (~75 Mo)


class Fenetre(QMainWindow):
    def __init__(self):
        super().__init__()
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
        btn_ncbi = QPushButton("Télécharger")
        btn_ncbi.clicked.connect(self.telecharger)
        for w in (self.accession, self.email, self.cible, btn_ncbi):
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
        self.mode = QComboBox()
        self.mode.addItems(["global", "local"])
        options.addWidget(self.mode)
        self.btn_align = QPushButton("Aligner")
        self.btn_align.setEnabled(False)
        self.btn_align.clicked.connect(self.lancer_alignement)
        options.addWidget(self.btn_align)
        options.addStretch(1)
        layout.addLayout(options)

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
        self.labels[cle].setText(f"{ident} — {len(seq)} pb{suite}")
        pret = all(self.seqs.values())
        self.btn_dot.setEnabled(pret)
        self.btn_align.setEnabled(pret)

    def tracer(self):
        (id_a, a), (id_b, b) = self.seqs["A"], self.seqs["B"]
        if len(a) * len(b) > MAX_CASES:
            QMessageBox.warning(self, "Trop grand",
                                "Séquences trop longues pour un dotplot complet.")
            return
        try:
            matrice = dotplot(a, b, fenetre=self.spin.value())
        except ValueError as e:
            QMessageBox.warning(self, "Erreur", str(e))
            return
        self.figure.clear()
        tracer_dotplot(matrice, self.figure.add_subplot(111), id_a, id_b)
        self.canvas.draw()

    def lancer_alignement(self):
        (id_a, a), (id_b, b) = self.seqs["A"], self.seqs["B"]
        if len(a) * len(b) > MAX_CASES:
            QMessageBox.warning(self, "Trop grand",
                                "Séquences trop longues pour cet alignement.")
            return
        mode = self.mode.currentText()
        score, la, lb = aligner(a.upper(), b.upper(), mode=mode)
        milieu = "".join("|" if x == y else " " for x, y in zip(la, lb))
        ident = milieu.count("|")
        texte = [f"Mode : {mode} | score : {score} | identité : {ident}/{len(la)}", ""]
        for i in range(0, len(la), 60):
            texte += [la[i:i + 60], milieu[i:i + 60], lb[i:i + 60], ""]
        self.sortie.setPlainText("\n".join(texte))

    def telecharger(self):
        email = self.email.text().strip()
        QSettings("bioapp", "bioapp").setValue("email", email)
        try:
            chemin = telecharger_fasta(self.accession.text(), email)
        except Exception as e:
            QMessageBox.warning(self, "NCBI", f"Téléchargement impossible : {e}")
            return
        self.charger_chemin(self.cible.currentText(), str(chemin))
