from gui.couleurs import blocs_alignement, colorer_alignement, milieu_alignement


def test_milieu():
    assert milieu_alignement("AC-T", "AGGT") == "|  |"


def test_blocs_largeur():
    assert blocs_alignement("ACGT", "ACGT", largeur=2) == [
        "AC", "||", "AC", "", "GT", "||", "GT", ""]


def test_couleurs_identite_difference_gap():
    html = colorer_alignement(["entête", ""] + blocs_alignement("AC-T", "AGGT"))
    assert "entête" in html
    assert '<span style="color:#2e9e4f">A</span>' in html   # identité
    assert '<span style="color:#d6453d">C</span>' in html   # différence
    assert '<span style="color:#d98a00">-</span>' in html   # gap


def test_html_echappe():
    assert "&lt;" in colorer_alignement(["a<b", ""])
