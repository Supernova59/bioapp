import pytest
from core.comparaison import comparer, pourcentage_identite, choisir_bareme


def test_trois_sequences_donnent_trois_paires():
    r = comparer({"a": "ACGTACGTAC", "b": "ACGTACGTAC", "c": "ACGAACGTTC"})
    assert [(x["a"], x["b"]) for x in r] == [("a", "b"), ("a", "c"), ("b", "c")]
    assert r[0]["score"] == 10.0 and r[0]["identite"] == 100.0
    assert r[1]["score"] == 6.0 and r[1]["identite"] == 80.0


def test_deux_sequences_donnent_une_paire():
    assert len(comparer({"x": "ACGT", "y": "ACGT"})) == 1


def test_nombre_de_sequences_invalide():
    with pytest.raises(ValueError):
        comparer({"x": "ACGT"})


def test_proteines_utilisent_blosum62():
    r = comparer({"p1": "HEAGAWGHEE", "p2": "PAWHEAE"}, mode="local")[0]
    assert r["score"] == 18.0 and "BLOSUM62" in r["bareme"]


def test_identite_ignore_les_gaps():
    assert pourcentage_identite("AC-T", "ACGT") == 75.0
    assert pourcentage_identite("", "") == 0.0


def test_bareme_adn_modifiable():
    assert choisir_bareme("ACGT", "ACGT", match=2)[0]["match"] == 2
