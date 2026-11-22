import time

from PySide6.QtCore import QTimer
from PySide6.QtGui import QFontDatabase
from PySide6.QtWidgets import (QAbstractItemView, QApplication, QComboBox, QDialog,
                               QHBoxLayout, QLabel, QMessageBox, QPlainTextEdit,
                               QPushButton, QSpinBox, QTableWidget, QTableWidgetItem,
                               QVBoxLayout)

from core.blast import PROGRAMMES, lancer_blast
from gui.travailleur import Travailleur

EVALUES = ["10", "1", "0.01", "1e-5", "1e-10"]


def _texte_alignement(h, largeur=60):
    lignes = [f"{h['accession']} — {h['titre']}",
              f"Requête {h['debut_requete']}–{h['fin_requete']} | "
              f"cible {h['debut_cible']}–{h['fin_cible']} | score {h['score']} "
              f"({h['bits']} bits) | identité {h['identite']:.1f} %", ""]
    for i in range(0, len(h["requete"]), largeur):
        lignes += [h["requete"][i:i + largeur], h["milieu"][i:i + largeur],
                   h["cible"][i:i + largeur], ""]
    return "\n".join(lignes)


class DialogueBlast(QDialog):
    def __init__(self, fenetre):
        super().__init__(fenetre)
        self.fenetre = fenetre
        self.tache = None
        self.hits = []
        self.debut = 0.0
        self.setWindowTitle("BLAST (serveur NCBI)")
        self.resize(760, 600)
        layout = QVBoxLayout(self)

        ligne = QHBoxLayout()
        self.cible = QComboBox()
        self.cible.addItems(["A", "B"])
        self.programme = QComboBox()
        self.programme.addItems(list(PROGRAMMES))
        self.programme.setCurrentText("blastp")
        self.base = QComboBox()
        self.programme.currentTextChanged.connect(self._bases)
        self._bases("blastp")
        self.evalue = QComboBox()
        self.evalue.addItems(EVALUES)
        self.evalue.setCurrentText("0.01")
        self.maximum = QSpinBox()
        self.maximum.setRange(1, 50)
        self.maximum.setValue(10)
        self.btn = QPushButton("Lancer BLAST")
        self.btn.clicked.connect(self.lancer)
        for w in (QLabel("Séquence :"), self.cible, self.programme, self.base,
                  QLabel("E ≤"), self.evalue, QLabel("hits"), self.maximum, self.btn):
            ligne.addWidget(w)
        layout.addLayout(ligne)

        self.etat = QLabel("blastn/blastx : séquence d'ADN ; blastp : protéine. "
                           "Une requête peut durer plusieurs minutes.")
        layout.addWidget(self.etat)
        self.table = QTableWidget(0, 6)
        self.table.setHorizontalHeaderLabels(
            ["Accession", "Description", "E-value", "Identité %", "Couverture %", "Score"])
        self.table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.itemSelectionChanged.connect(self.montrer)
        layout.addWidget(self.table, 1)
        self.detail = QPlainTextEdit()
        self.detail.setReadOnly(True)
        self.detail.setFont(QFontDatabase.systemFont(QFontDatabase.FixedFont))
        self.detail.setMaximumHeight(190)
        layout.addWidget(self.detail)

        self.minuteur = QTimer(self)
        self.minuteur.timeout.connect(self._tic)
        QApplication.instance().aboutToQuit.connect(self._arreter)

    def _bases(self, programme):
        self.base.clear()
        self.base.addItems(list(PROGRAMMES.get(programme, ())))

    def _tic(self):
        self.etat.setText(f"Requête envoyée à NCBI… {int(time.time() - self.debut)} s "
                          "(plusieurs minutes possibles)")

    def lancer(self):
        if self.tache is not None and self.tache.isRunning():
            return
        entree = self.fenetre.seqs.get(self.cible.currentText())
        if entree is None:
            QMessageBox.information(self, "BLAST", "Charge d'abord la séquence "
                                    f"{self.cible.currentText()}.")
            return
        self.btn.setEnabled(False)
        self.table.setRowCount(0)
        self.detail.clear()
        self.debut = time.time()
        self.minuteur.start(1000)
        self._tic()
        self.tache = Travailleur(lancer_blast, entree[1], self.programme.currentText(),
                                 self.base.currentText(), float(self.evalue.currentText()),
                                 self.maximum.value())
        self.tache.termine.connect(self._fin)
        self.tache.echec.connect(self._echec)
        self.tache.start()

    def _termine(self):
        self.minuteur.stop()
        self.tache.wait()
        self.btn.setEnabled(True)

    def _echec(self, message):
        self._termine()
        self.etat.setText(f"Erreur : {message}")

    def _fin(self, hits):
        self._termine()
        self.hits = hits
        self.etat.setText(f"{len(hits)} hit(s) en {int(time.time() - self.debut)} s"
                          if hits else "Aucun hit sous ce seuil.")
        self.table.setRowCount(len(hits))
        for i, h in enumerate(hits):
            valeurs = [h["accession"], h["titre"], f"{h['evalue']:.2g}",
                       f"{h['identite']:.1f}", f"{h['couverture']:.1f}", f"{h['score']}"]
            for j, v in enumerate(valeurs):
                self.table.setItem(i, j, QTableWidgetItem(v))
        self.table.resizeColumnsToContents()
        if hits:
            self.table.selectRow(0)

    def montrer(self):
        lignes = self.table.selectionModel().selectedRows()
        if lignes and lignes[0].row() < len(self.hits):
            self.detail.setPlainText(_texte_alignement(self.hits[lignes[0].row()]))

    def _arreter(self):
        if self.tache is not None and self.tache.isRunning():
            self.tache.terminate()
            self.tache.wait(1000)

    def done(self, resultat):
        if self.tache is not None and self.tache.isRunning():
            self.hide()          # la requête continue ; rouvrir le bouton BLAST la retrouve
            return
        super().done(resultat)


class BlastMixin:
    def ouvrir_blast(self):
        if getattr(self, "_blast", None) is None:
            self._blast = DialogueBlast(self)
        self._blast.show()
        self._blast.raise_()
