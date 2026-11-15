import numpy as np
import pytest
from core.dotplot import dotplot


def test_identique_donne_diagonale():
    assert np.array_equal(dotplot("ACGT", "ACGT"), np.eye(4, dtype=bool))


def test_forme_avec_fenetre():
    assert dotplot("ACGTAC", "ACGT", fenetre=3).shape == (4, 2)


def test_fenetre_filtre_le_bruit():
    m = dotplot("ACGT", "ACGT", fenetre=2)
    assert np.array_equal(m, np.eye(3, dtype=bool))


def test_fenetre_trop_grande():
    with pytest.raises(ValueError):
        dotplot("AC", "ACGT", fenetre=3)
