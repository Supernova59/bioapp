"""Cas tirés des TP du cours « Sequence comparison » (valeurs vérifiées à la main ou par Biopython)."""
from core.alignement import aligner_complet
from core.comparaison import comparer
from core.sequences import complement_inverse
from core.significativite import significativite

BAREME_COURS = dict(match=5, mismatch=-4, ouverture=-10, extension=-0.5)

DNA_H = ("atgagtcactatattgaattaactgaagaaaattttgaaagcaccattaaaaaaggggttgcgttagtggatttttgggcaccatggtgtgg"
         "tccttgtaagatgctatcccctgtgattgatgaattagctagcgaatatgaaggtaaggctaagatttgtaaagttaataccgatgagcaag"
         "aagaattgagcgcgaaatttggtattaggagcattcctacgcttttattcacaaaagatggcgaagttgtccatcagttggtgggcgtgcaaa"
         "ctaaagtcgctttaaaagagcaattgaacaagcttttaggctag")
DNA_S = ("atggcaatcgtaaaagtaacagatgcagattttgattcaaaagtagaatctggtgtacaactagtagatttttgggcaacatggtgtggtccatg"
         "taaaatgatcgctccggtattagaagaattagcagctgactatgaaggtaaagctgacattttaaaattagatgttgatgaaaatccatcaact"
         "gcagctaaatatgaagtgatgagtattccaacattaatcgtctttaaagacggtcaaccagttgataaagttgttggtttccaaccaaaagaaaa"
         "cttagctgaagttttagataaacatttataa")
PROT_H = "MSHYIELTEENFESTIKKGVALVDFWAPWCGPCKMLSPVIDELASEYEGKAKICKVNTDEQEELSAKFGIRSIPTLLFTKDGEVVHQLVGVQTKVALKEQLNKLLG"
PROT_S = "MAIVKVTDADFDSKVESGVQLVDFWATWCGPCKMIAPVLEELAADYEGKADILKLDVDENPSTAAKYEVMSIPTLIVFKDGQPVDKVVGFQPKENLAEVLDKHL"


def test_alignement_manuel_semi_global():
    r = aligner_complet("AAGTCATTCGCACATCG", "AACACATCG", "semi-global", **BAREME_COURS)
    assert r[:3] == (31.5, "AAGTCATTCGCACATCG", "AA--------CACATCG")


def test_aagg_ag_semi_global_contre_global():
    p = dict(match=1, mismatch=-1, ouverture=-10, extension=-5)
    assert aligner_complet("AAGG", "AG", "semi-global", **p)[0] == 2.0
    assert aligner_complet("AAGG", "AG", "global", **p)[0] == -13.0


def test_thioredoxine_adn_local():
    r = comparer({"H": DNA_H, "S": DNA_S}, mode="local", **BAREME_COURS)[0]
    assert r["score"] == 556.5 and round(r["identite"], 1) == 60.8


def test_thioredoxine_proteine_local():
    r = comparer({"H": PROT_H, "S": PROT_S}, mode="local")[0]
    assert r["score"] == 290.0 and round(r["identite"], 1) == 51.5


def test_significativite_proteine_plus_forte_que_adn():
    adn = significativite(DNA_H, DNA_S, melanges=200, graine=1, **BAREME_COURS)
    prot = significativite(PROT_H, PROT_S, melanges=200, graine=1,
                           ouverture=-10, extension=-1, matrice="BLOSUM62")
    assert adn["meilleurs"] == 0 and prot["meilleurs"] == 0
    assert prot["z"] > adn["z"] > 5


def test_significativite_sequences_sans_rapport():
    r = significativite("ACGTACGTAC", "TTTTTTTTTT", melanges=50, graine=1)
    assert r["p"] == 1.0


def test_complement_inverse():
    assert complement_inverse("ATGCN") == "NGCAT"
    s = "tcgcgcgcgcgctgtaagagtgtgtgtggct"
    assert complement_inverse(complement_inverse(s)) == s
