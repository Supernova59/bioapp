import sys
from PySide6.QtWidgets import QApplication
from gui.fenetre import Fenetre

app = QApplication(sys.argv)
fenetre = Fenetre()
fenetre.resize(700, 800)
fenetre.show()
sys.exit(app.exec())
