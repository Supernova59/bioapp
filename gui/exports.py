from pathlib import Path
from PySide6.QtWidgets import QFileDialog, QMessageBox
from core.export import fasta_aligne


class ExportMixin:
    def _chemin_sortie(self, titre, filtre, extension):
        chemin, _ = QFileDialog.getSaveFileName(self, titre, "", filtre)
        if not chemin:
            return None
        if "." not in chemin.rsplit("/", 1)[-1]:
            chemin += extension
        return chemin

    def exporter_dotplot(self):
        if not self.figure.axes:
            QMessageBox.information(self, "Export", "Aucun dotplot à exporter.")
            return
        chemin = self._chemin_sortie("Exporter le dotplot",
                                     "PNG (*.png);;SVG (*.svg)", ".png")
        if chemin:
            self.figure.savefig(chemin, dpi=200, bbox_inches="tight")
            self.statusBar().showMessage(f"Dotplot exporté : {chemin}", 5000)

    def exporter_texte(self):
        texte = self.sortie.toPlainText()
        if not texte:
            QMessageBox.information(self, "Export", "Aucun alignement à exporter.")
            return
        chemin = self._chemin_sortie("Exporter l'alignement", "Texte (*.txt)", ".txt")
        if chemin:
            Path(chemin).write_text(texte + "\n", encoding="utf-8")
            self.statusBar().showMessage(f"Alignement exporté : {chemin}", 5000)

    def exporter_fasta(self):
        dernier = getattr(self, "dernier", None)
        if not dernier or not dernier[1]:
            QMessageBox.information(self, "Export", "Aucun alignement à exporter.")
            return
        (id_a, id_b), la, lb = dernier
        chemin = self._chemin_sortie("Exporter en FASTA aligné",
                                     "FASTA (*.fasta *.fa)", ".fasta")
        if chemin:
            Path(chemin).write_text(fasta_aligne(id_a, la, id_b, lb), encoding="utf-8")
            self.statusBar().showMessage(f"FASTA exporté : {chemin}", 5000)
