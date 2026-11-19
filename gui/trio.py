from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure
from PySide6.QtCore import QSettings
from PySide6.QtGui import QColor, QFontDatabase
from PySide6.QtWidgets import (QAbstractItemView, QComboBox, QDialog, QFileDialog,
                               QHBoxLayout, QLabel, QMessageBox, QPushButton,
                               QTableWidget, QTableWidgetItem, QTabWidget, QTextEdit,
                               QVBoxLayout, QWidget)

from core.comparaison import comparer
from core.dotplot import dotplot, tracer_dotplot
from core.fasta import lire_fasta
from gui.couleurs import blocs_alignement, colorer_alignement
from gui.travailleur import Travailleur

MAX_CASES = 25_000_000


def _calculer(sequences, mode, bareme, fenetre_dot):
    resultats = comparer(sequences, mode=mode, **bareme)
    matrices = {}
    for r in resultats:
        a, b = sequences[r["a"]].upper(), sequences[r["b"]].upper()
        if len(a) * len(b) <= MAX_CASES:
            try:
                matrices[(r["a"], r["b"])] = dotplot(a, b, fenetre=fenetre_dot)
            except ValueError:
                pass
    return resultats, matrices


class DialogueTrio(QDialog):
    def __init__(self, fenetre):
        super().__init__(fenetre)
        self.fenetre = fenetre
        self.tache = None
        self.seq_c = None
        self.setWindowTitle("Comparaison de 3 séquences (A, B, C)")
        self.resize(760, 640)
        layout = QVBoxLayout(self)

        ligne = QHBoxLayout()
        btn_c = QPushButton("Charger C")
        btn_c.clicked.connect(self.charger_c)
        self.label_c = QLabel("C : (rien chargé) — A et B viennent de la fenêtre principale")
        self.mode = QComboBox()
        self.mode.addItems(["global", "semi-global", "local"])
        self.btn_go = QPushButton("Comparer")
        self.btn_go.clicked.connect(self.comparer)
        for w in (btn_c, self.label_c):
            ligne.addWidget(w)
        ligne.addStretch(1)
        ligne.addWidget(self.mode)
        ligne.addWidget(self.btn_go)
        layout.addLayout(ligne)

        self.onglets = QTabWidget()
        layout.addWidget(self.onglets, 1)

    def charger_c(self):
        reglages = QSettings("bioapp", "bioapp")
        chemin, _ = QFileDialog.getOpenFileName(
            self, "Ouvrir le FASTA C", reglages.value("dossier", ""),
            "Tous (*);;FASTA (*.fasta *.fa *.fna *.faa)")
        if not chemin:
            return
        try:
            entrees = lire_fasta(chemin)
            if not entrees:
                raise ValueError("aucune séquence trouvée")
        except Exception as e:
            QMessageBox.warning(self, "Erreur", f"Lecture impossible : {e}")
            return
        self.seq_c = entrees[0]
        self.label_c.setText(f"C : {entrees[0][0]} — {len(entrees[0][1])}")

    def comparer(self):
        if self.tache is not None and self.tache.isRunning():
            return
        seqs = self.fenetre.seqs
        if not (seqs.get("A") and seqs.get("B") and self.seq_c):
            QMessageBox.information(self, "Comparaison",
                                    "Charge A et B dans la fenêtre principale, puis C ici.")
            return
        sequences = {f"A ({seqs['A'][0]})": seqs["A"][1],
                     f"B ({seqs['B'][0]})": seqs["B"][1],
                     f"C ({self.seq_c[0]})": self.seq_c[1]}
        bareme = {k: v.value() for k, v in self.fenetre.params.items()}
        self.btn_go.setEnabled(False)
        self.onglets.clear()
        self.onglets.addTab(QLabel("Calcul en cours…"), "Résumé")
        self.tache = Travailleur(_calculer, sequences, self.mode.currentText(),
                                 bareme, self.fenetre.spin.value())
        self.tache.termine.connect(self._fin)
        self.tache.echec.connect(self._echec)
        self.tache.start()

    def _echec(self, message):
        self.tache.wait()
        self.btn_go.setEnabled(True)
        self.onglets.clear()
        self.onglets.addTab(QLabel(f"Erreur : {message}"), "Résumé")

    def _fin(self, resultat):
        self.tache.wait()
        self.btn_go.setEnabled(True)
        resultats, matrices = resultat
        self.onglets.clear()
        self.onglets.addTab(self._page_resume(resultats), "Résumé")
        for r in resultats:
            self.onglets.addTab(self._page_paire(r, matrices.get((r["a"], r["b"]))),
                                f"{r['a'][0]}–{r['b'][0]}")

    def _page_resume(self, resultats):
        page = QWidget()
        layout = QVBoxLayout(page)
        meilleur = max(resultats, key=lambda r: (r["identite"], r["score"]))
        table = QTableWidget(len(resultats), 4)
        table.setHorizontalHeaderLabels(["Paire", "Score", "Identité (%)", "Colonnes"])
        table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        table.horizontalHeader().setStretchLastSection(True)
        for i, r in enumerate(resultats):
            valeurs = [f"{r['a']}  vs  {r['b']}", f"{r['score']}",
                       f"{r['identite']:.1f}", f"{len(r['ligne_a'])}"]
            for j, v in enumerate(valeurs):
                item = QTableWidgetItem(v)
                if r is meilleur:
                    item.setBackground(QColor("#2e9e4f"))
                    item.setForeground(QColor("#ffffff"))
                table.setItem(i, j, item)
        table.resizeColumnsToContents()
        layout.addWidget(table)
        layout.addWidget(QLabel(
            f"Paire la plus proche : {meilleur['a']} vs {meilleur['b']} "
            f"({meilleur['identite']:.1f} % d'identité)\n{resultats[0]['bareme']}"))
        layout.addStretch(1)
        return page

    def _page_paire(self, r, matrice):
        page = QWidget()
        layout = QVBoxLayout(page)
        if matrice is not None:
            figure = Figure(figsize=(4, 4))
            tracer_dotplot(matrice, figure.add_subplot(111), r["a"], r["b"])
            layout.addWidget(FigureCanvasQTAgg(figure), 1)
        else:
            layout.addWidget(QLabel("Dotplot non calculé (séquences trop longues "
                                    "ou fenêtre plus grande qu'une séquence)."))
        zone = QTextEdit()
        zone.setReadOnly(True)
        zone.setFont(QFontDatabase.systemFont(QFontDatabase.FixedFont))
        zone.setMaximumHeight(200)
        entete = [f"{r['a']} vs {r['b']} | score : {r['score']} | "
                  f"identité : {r['identite']:.1f} %", ""]
        zone.setHtml(colorer_alignement(entete + blocs_alignement(r["ligne_a"], r["ligne_b"])))
        layout.addWidget(zone)
        return page

    def done(self, resultat):
        if self.tache is not None and self.tache.isRunning():
            self.tache.wait(10000)
        super().done(resultat)


class TrioMixin:
    def ouvrir_trio(self):
        self._trio = DialogueTrio(self)
        self._trio.show()
