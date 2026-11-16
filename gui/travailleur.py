from PySide6.QtCore import QThread, Signal


class Travailleur(QThread):
    """Exécute une fonction en arrière-plan et renvoie son résultat ou l'erreur."""
    termine = Signal(object)
    echec = Signal(str)

    def __init__(self, fonction, *args):
        super().__init__()
        self.fonction = fonction
        self.args = args

    def run(self):
        try:
            self.termine.emit(self.fonction(*self.args))
        except Exception as e:
            self.echec.emit(str(e))
