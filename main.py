import sys
from PySide6.QtWidgets import QApplication
from gui.fenetre import Fenetre
from gui.presence import Presence

app = QApplication(sys.argv)
presence = Presence()
presence.maj("Analyse de séquences")
app.aboutToQuit.connect(presence.fermer)
fenetre = Fenetre(presence)
fenetre.resize(700, 800)
fenetre.show()
sys.exit(app.exec())
