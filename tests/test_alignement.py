import pytest
from core.alignement import aligner


def test_global_identique():
    assert aligner("ACGT", "ACGT") == (4, "ACGT", "ACGT")


def test_global_gap():
    assert aligner("ACGT", "AGT") == (1, "ACGT", "A-GT")


def test_local():
    score, a, b = aligner("TTACGTT", "GGACGGG", mode="local")
    assert (score, a, b) == (3, "ACG", "ACG")


def test_mode_invalide():
    with pytest.raises(ValueError):
        aligner("A", "A", mode="xyz")
