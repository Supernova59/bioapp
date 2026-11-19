import pytest
from core.annotation import chercher_orfs, profil_hydropathie, regions_hydrophobes, traduire
from core.sequences import complement_inverse


def test_hydropathie_valeurs_kd():
    assert profil_hydropathie("AI", fenetre=1) == [(1, 1.8), (2, 4.5)]


def test_hydropathie_moyenne_et_inconnu():
    assert profil_hydropathie("AIA", fenetre=3) == [(2, (1.8 + 4.5 + 1.8) / 3)]
    with pytest.raises(ValueError):
        profil_hydropathie("AZA", fenetre=1)


def test_hydropathie_fenetre_invalide():
    with pytest.raises(ValueError):
        profil_hydropathie("AAAA", fenetre=2)
    with pytest.raises(ValueError):
        profil_hydropathie("AA", fenetre=9)


def test_region_transmembranaire():
    prot = "D" * 15 + "L" * 25 + "D" * 15
    r = regions_hydrophobes(prot)
    assert len(r) == 1 and r[0][0] >= 11 and r[0][1] <= 45
    assert regions_hydrophobes("D" * 40) == []


def test_traduire_cadres():
    assert traduire("ATGAAATTTGGGCCCTAA") == "MKFGP*"
    assert traduire("CATGAAATTTGGGCCCTAA", cadre=2) == "MKFGP*"
    assert traduire("ATGAAATTTGGGCCCTAA"[:-1]) == "MKFGP"   # codon final incomplet ignoré
    assert traduire(complement_inverse("ATGAAATTTGGGCCCTAA"), cadre=-1) == "MKFGP*"


def test_orf_simple_et_decalage():
    orf = "ATGAAATTTGGGCCCTAA"
    r = chercher_orfs("CC" + orf + "GG", longueur_min=3)
    assert [(o["brin"], o["debut"], o["fin"], o["proteine"]) for o in r] == [("+", 3, 20, "MKFGP")]


def test_orf_brin_inverse():
    orf = "ATGAAATTTGGGCCCTAA"
    adn = complement_inverse("CC" + orf + "GG")
    r = chercher_orfs(adn, longueur_min=3)
    assert [(o["brin"], o["debut"], o["fin"]) for o in r] == [("-", 3, 20)]
    assert chercher_orfs(adn, longueur_min=3, deux_brins=False) == []


def test_orf_sans_stop_ignore_et_longueur_min():
    assert chercher_orfs("ATGAAATTTGGGCCC", longueur_min=3) == []
    assert chercher_orfs("ATGAAATTTGGGCCCTAA", longueur_min=30) == []
