from PySide6.QtGui import QFontDatabase
from PySide6.QtWidgets import (QCheckBox, QComboBox, QDialog, QHBoxLayout, QLabel,
                               QLineEdit, QMessageBox, QPlainTextEdit, QPushButton,
                               QSpinBox, QTabWidget, QVBoxLayout, QWidget)

from core.comparaison import choisir_bareme
from core.motifs import chercher_motif
from core.significativite import significativite
from gui.travailleur import Travailleur


def _zone_texte():
    zone = QPlainTextEdit()
    zone.setReadOnly(True)
    zone.setFont(QFontDatabase.systemFont(QFontDatabase.FixedFont))
    return zone


class DialogueOutils(QDialog):
    def __init__(self, fenetre):
        super().__init__(fenetre)
        self.fenetre = fenetre
        self.tache = None
        self.setWindowTitle("Outils")
        self.resize(600, 460)
        onglets = QTabWidget()
        onglets.addTab(self._onglet_motifs(), "Motifs (Prosite)")
        onglets.addTab(self._onglet_signif(), "Significativité")
        from gui.outils2 import page_hydro, page_msa, page_orf
        onglets.addTab(page_orf(self), "Traduction / ORF")
        onglets.addTab(page_hydro(self), "Hydropathie")
        onglets.addTab(page_msa(self), "Alignement multiple")
        from gui.outils2 import page_hydro, page_msa, page_orf
        onglets.addTab(page_orf(self), "Traduction / ORF")
        onglets.addTab(page_hydro(self), "Hydropathie")
        onglets.addTab(page_msa(self), "Alignement multiple")
        QVBoxLayout(self).addWidget(onglets)

    # ---------- utilitaires ----------
    def _sequence(self, cle):
        entree = self.fenetre.seqs.get(cle)
        if entree is None:
            QMessageBox.information(self, "Outils", f"Charge d'abord la séquence {cle}.")
        return entree

    # ---------- motifs ----------
    def _onglet_motifs(self):
        page = QWidget()
        layout = QVBoxLayout(page)
        ligne = QHBoxLayout()
        self.motif = QLineEdit()
        self.motif.setPlaceholderText("motif Prosite, ex. C-x(2,4)-C-x(3)-[LIVMFYWC]-x(8)-H-x(3,5)-H")
        self.cible_motif = QComboBox()
        self.cible_motif.addItems(["A", "B"])
        self.deux_brins = QCheckBox("Deux brins (ADN)")
        bouton = QPushButton("Chercher")
        bouton.clicked.connect(self.chercher)
        self.motif.returnPressed.connect(self.chercher)
        for w in (self.motif, self.cible_motif, self.deux_brins, bouton):
            ligne.addWidget(w)
        layout.addLayout(ligne)
        self.res_motifs = _zone_texte()
        layout.addWidget(self.res_motifs)
        return page

    def chercher(self):
        cle = self.cible_motif.currentText()
        entree = self._sequence(cle)
        if entree is None:
            return
        ident, seq = entree
        try:
            trouves = chercher_motif(seq, self.motif.text(), self.deux_brins.isChecked())
        except ValueError as e:
            QMessageBox.warning(self, "Motif", str(e))
            return
        lignes = [f"{ident} : {len(trouves)} occurrence(s)", ""]
        lignes += [f"{t['brin']}  {t['debut']}–{t['fin']}  {t['sequence']}" for t in trouves]
        self.res_motifs.setPlainText("\n".join(lignes))

    # ---------- significativité ----------
    def _onglet_signif(self):
        page = QWidget()
        layout = QVBoxLayout(page)
        ligne = QHBoxLayout()
        ligne.addWidget(QLabel("Mélanges :"))
        self.nb = QSpinBox()
        self.nb.setRange(10, 5000)
        self.nb.setValue(500)
        ligne.addWidget(self.nb)
        self.mode_signif = QComboBox()
        self.mode_signif.addItems(["local", "global", "semi-global"])
        ligne.addWidget(self.mode_signif)
        self.btn_calc = QPushButton("Calculer")
        self.btn_calc.clicked.connect(self.calculer)
        ligne.addWidget(self.btn_calc)
        ligne.addStretch(1)
        layout.addLayout(ligne)
        self.res_signif = _zone_texte()
        layout.addWidget(self.res_signif)
        return page

    @staticmethod
    def _calcul(a, b, n, mode, params):
        return significativite(a, b, melanges=n, mode=mode, **params)

    def calculer(self):
        if self.tache is not None and self.tache.isRunning():
            return
        ea, eb = self._sequence("A"), self._sequence("B")
        if ea is None or eb is None:
            return
        (id_a, a), (id_b, b) = ea, eb
        p = {k: v.value() for k, v in self.fenetre.params.items()}
        params, description = choisir_bareme(a.upper(), b.upper(), **p)
        self.contexte = (id_a, id_b, self.mode_signif.currentText(), description)
        self.btn_calc.setEnabled(False)
        self.res_signif.setPlainText("Calcul en cours…")
        self.tache = Travailleur(self._calcul, a, b, self.nb.value(),
                                 self.mode_signif.currentText(), params)
        self.tache.termine.connect(self._fin_signif)
        self.tache.echec.connect(self._echec_signif)
        self.tache.start()

    def _fin_signif(self, d):
        self.tache.wait()
        self.btn_calc.setEnabled(True)
        id_a, id_b, mode, description = self.contexte
        lignes = [
            f"{id_a} vs {id_b} | mode {mode} | {description}",
            "",
            f"Score réel : {d['score']}",
            f"Scores mélangés ({d['melanges']}) : moyenne {d['moyenne']:.1f}, "
            f"écart-type {d['ecart_type']:.1f}",
            f"z-score : {d['z']:.1f}",
            f"Mélanges au moins aussi bons : {d['meilleurs']}/{d['melanges']}",
            f"p empirique : {d['p']:.4f}",
        ]
        if d["meilleurs"] == 0:
            lignes.append(f"(p ne peut pas descendre sous 1/{d['melanges'] + 1} : "
                          "augmente les mélanges pour affiner)")
        self.res_signif.setPlainText("\n".join(lignes))

    def _echec_signif(self, message):
        self.tache.wait()
        self.btn_calc.setEnabled(True)
        self.res_signif.setPlainText(f"Erreur : {message}")

    def done(self, resultat):
        if self.tache is not None and self.tache.isRunning():
            self.tache.wait(10000)
        page = getattr(self, "_page_msa", None)
        if page is not None:
            page.attendre()
        page = getattr(self, "_page_msa", None)
        if page is not None:
            page.attendre()
        super().done(resultat)


class OutilsMixin:
    def ouvrir_outils(self):
        self._outils = DialogueOutils(self)
        self._outils.show()
