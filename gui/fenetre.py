from PySide6.QtWidgets import (QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
                               QPushButton, QLabel, QFileDialog, QMessageBox, QSpinBox)
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure

from core.fasta import lire_fasta
from core.dotplot import dotplot, tracer_dotplot

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
        options.addStretch(1)
        layout.addLayout(options)

        self.figure = Figure(figsize=(5, 5))
        self.canvas = FigureCanvasQTAgg(self.figure)
        layout.addWidget(self.canvas, 1)

    def charger(self, cle):
        chemin, _ = QFileDialog.getOpenFileName(
            self, "Ouvrir un FASTA", "", "FASTA (*.fasta *.fa *.fna *.faa);;Tous (*)")
        if not chemin:
            return
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
        self.btn_dot.setEnabled(all(self.seqs.values()))

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
