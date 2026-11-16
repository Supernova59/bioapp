import pytest
from core.motifs import chercher_motif, prosite_vers_regex


def test_conversion():
    assert prosite_vers_regex("C-x(2,4)-[LIV]-{P}-G.") == "C.{2,4}[LIV][^P]G"


def test_doigt_de_zinc_c2h2():
    motif = "C-x(2,4)-C-x(3)-[LIVMFYWC]-x(8)-H-x(3,5)-H"
    seq = "MKK" + "CAAC" + "AAA" + "L" + "A" * 8 + "H" + "AAA" + "H" + "GG"
    r = chercher_motif(seq, motif)
    assert len(r) == 1 and r[0]["debut"] == 4 and r[0]["fin"] == 24


def test_ancres_debut_et_fin():
    assert [r["debut"] for r in chercher_motif("MASMAS", "<M-A-S")] == [1]
    assert [r["debut"] for r in chercher_motif("GAGAG", "G>")] == [5]
    assert [r["debut"] for r in chercher_motif("AGA", "A-[G>]")] == [1, 3]


def test_occurrences_chevauchantes():
    assert [r["debut"] for r in chercher_motif("AAAA", "A-A")] == [1, 2, 3]


def test_deux_brins():
    r = chercher_motif("AAAGGGTTT", "A-A-A", deux_brins=True)
    assert [(x["debut"], x["fin"], x["brin"]) for x in r] == [(1, 3, "+"), (7, 9, "-")]


def test_motif_invalide():
    for mauvais in ("", "C-x(2,", "1-2", "C-[ab]"):
        with pytest.raises(ValueError):
            chercher_motif("ACGT", mauvais)
