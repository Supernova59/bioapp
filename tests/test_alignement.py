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


def test_local_bat_global_sur_bords_differents():
    a = "T" * 10 + "ATGCGTACGTTAGC" + "T" * 10
    b = "G" * 10 + "ATGCGTACGTTAGC" + "C" * 9
    assert aligner(a, b, mode="local")[0] == 14
    assert aligner(a, b, mode="global")[0] == -7


from core.alignement import aligner_complet


def test_positions_locales():
    a = "T" * 10 + "ATGCGTACGTTAGC" + "T" * 10
    b = "G" * 10 + "ATGCGTACGTTAGC" + "C" * 9
    score, _, _, pos_a, pos_b = aligner_complet(a, b, mode="local")
    assert score == 14 and pos_a == (10, 24) and pos_b == (10, 24)


def test_bareme_modifiable():
    assert aligner_complet("ACGT", "ACGT", match=2)[0] == 8


def test_local_sans_ressemblance():
    r = aligner_complet("AAAA", "TTTT", mode="local")
    assert r[0] == 0 and r[1] == "" and r[2] == ""


from core.alignement import type_sequence


def test_type_sequence():
    assert type_sequence("ACGTNacgt") == "nucleotide"
    assert type_sequence("MKVLAAGIVG") == "proteine"


def test_blosum62_identique():
    r = aligner_complet("ACDEF", "ACDEF", matrice="BLOSUM62", ouverture=-10, extension=-1)
    assert r[0] == 30


def test_blosum62_local():
    r = aligner_complet("HEAGAWGHEE", "PAWHEAE", mode="local", matrice="BLOSUM62",
                        ouverture=-10, extension=-1)
    assert r[0] == 18


def test_semi_global_sequence_incluse():
    r = aligner_complet("GGGACGTACGTGGG", "ACGTACGT", mode="semi-global")
    assert r[:3] == (8.0, "GGGACGTACGTGGG", "---ACGTACGT---")
    assert aligner_complet("GGGACGTACGTGGG", "ACGTACGT")[0] == 0.0  # global : pénalisé


def test_semi_global_chevauchement():
    r = aligner_complet("TTTTACGTACGT", "ACGTACGT", mode="semi-global")
    assert r[:3] == (8.0, "TTTTACGTACGT", "----ACGTACGT")


def test_semi_global_identique():
    assert aligner_complet("ACGT", "ACGT", mode="semi-global")[0] == 4.0
