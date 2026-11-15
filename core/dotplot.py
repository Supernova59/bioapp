import numpy as np


def dotplot(seq_a, seq_b, fenetre=1, seuil=None):
    """Matrice booléenne (len(a)-fenetre+1, len(b)-fenetre+1).
    Une case est True si au moins `seuil` positions sur `fenetre` sont identiques
    le long de la diagonale (seuil = fenetre par défaut : correspondance exacte)."""
    if fenetre < 1:
        raise ValueError("fenetre doit être >= 1")
    seuil = fenetre if seuil is None else seuil
    a = np.frombuffer(seq_a.upper().encode(), dtype=np.uint8)
    b = np.frombuffer(seq_b.upper().encode(), dtype=np.uint8)
    n, m = len(a) - fenetre + 1, len(b) - fenetre + 1
    if n <= 0 or m <= 0:
        raise ValueError("fenetre plus grande qu'une des séquences")
    egal = a[:, None] == b[None, :]
    score = np.zeros((n, m), dtype=np.int16)
    for k in range(fenetre):
        score += egal[k:k + n, k:k + m]
    return score >= seuil


def tracer_dotplot(matrice, ax, nom_a="séquence A", nom_b="séquence B"):
    """Dessine la matrice sur un axe matplotlib (aucune dépendance à PySide6)."""
    ax.imshow(~matrice, cmap="gray", aspect="auto", interpolation="nearest", origin="lower")
    ax.set_xlabel(nom_b)
    ax.set_ylabel(nom_a)
