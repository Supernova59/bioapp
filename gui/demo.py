from pathlib import Path

from PySide6.QtWidgets import (QDialog, QHBoxLayout, QLabel, QListWidget, QMessageBox,
                               QPlainTextEdit, QPushButton, QVBoxLayout)

from gui.scenarios import SCENARIOS

DOSSIER_DEMO = Path.home() / ".cache" / "bioapp" / "demo"


class DialogueDemo(QDialog):
    """Charge un scénario de test dans la fenêtre principale et le lance."""

    def __init__(self, fenetre):
        super().__init__(fenetre)
        self.fenetre = fenetre
        self.setWindowTitle("Scénarios de test")
        self.resize(520, 420)
        layout = QVBoxLayout(self)
        layout.addWidget(QLabel("Choisis un scénario : il charge A et B, règle les "
                                "paramètres et lance le calcul dans la fenêtre principale."))
        self.liste = QListWidget()
        for sc in SCENARIOS:
            self.liste.addItem(sc["titre"])
        self.liste.currentRowChanged.connect(self.montrer)
        self.liste.itemDoubleClicked.connect(lambda _: self.lancer())
        layout.addWidget(self.liste, 1)
        layout.addWidget(QLabel("Ce que tu dois voir :"))
        self.attendu = QPlainTextEdit()
        self.attendu.setReadOnly(True)
        self.attendu.setMaximumHeight(110)
        layout.addWidget(self.attendu)
        bas = QHBoxLayout()
        bouton = QPushButton("Charger et lancer")
        bouton.clicked.connect(self.lancer)
        bas.addStretch(1)
        bas.addWidget(bouton)
        layout.addLayout(bas)
        self.liste.setCurrentRow(0)

    def montrer(self, ligne):
        if 0 <= ligne < len(SCENARIOS):
            self.attendu.setPlainText(SCENARIOS[ligne]["attendu"])

    def lancer(self):
        ligne = self.liste.currentRow()
        if ligne < 0:
            return
        f = self.fenetre
        if getattr(f, "tache", None) is not None and f.tache.isRunning():
            QMessageBox.information(self, "Scénarios", "Un calcul est déjà en cours.")
            return
        sc = SCENARIOS[ligne]
        DOSSIER_DEMO.mkdir(parents=True, exist_ok=True)
        for cle in ("A", "B"):
            ident, seq = sc[cle]
            chemin = DOSSIER_DEMO / f"{ident}.fasta"
            chemin.write_text(f">{ident}\n{seq}\n")
            f.charger_chemin(cle, str(chemin))
        manques = []
        for nom, valeur in sc["params"].items():
            try:
                f.params[nom].setValue(valeur)
            except TypeError:
                manques.append(nom)      # case entière : valeur décimale impossible
        f.mode.setCurrentText(sc["mode"])
        if "fenetre" in sc:
            f.spin.setValue(sc["fenetre"])
        if hasattr(f, "chk_rc"):
            f.chk_rc.setChecked(sc.get("rc", False))
        if manques:
            QMessageBox.warning(self, "Scénarios",
                                "Paramètre(s) non réglé(s) : " + ", ".join(manques))
        if sc["action"] == "dotplot":
            f.tracer()
        else:
            f.lancer_alignement()


class DemoMixin:
    def ouvrir_demo(self):
        self._demo = DialogueDemo(self)
        self._demo.show()
