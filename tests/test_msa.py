import stat
import pytest
from core import msa


def _faux_programme(tmp_path, nom):
    """Script qui imite clustalo : recopie l'entrée en ajoutant un gap à la 1re séquence."""
    script = tmp_path / nom
    script.write_text(
        "#!/usr/bin/env python3\n"
        "import sys\n"
        "a = sys.argv\n"
        "entree = a[a.index('-i') + 1]\n"
        "sortie = a[a.index('-o') + 1]\n"
        "lignes = open(entree).read().split('\\n')\n"
        "with open(sortie, 'w') as f:\n"
        "    for i, l in enumerate(lignes):\n"
        "        f.write((l + '-' if i == 1 else l) + '\\n')\n")
    script.chmod(script.stat().st_mode | stat.S_IEXEC)
    return str(script)


def test_aligner_multiple_avec_faux_programme(tmp_path):
    prog = _faux_programme(tmp_path, "clustalo")
    r = msa.aligner_multiple([("x y", "ACGT"), ("z", "ACGT")], programme=prog)
    assert r == [("x y", "ACGT-"), ("z", "ACGT")]


def test_programme_absent(monkeypatch):
    monkeypatch.setattr(msa, "trouver_programme", lambda: None)
    with pytest.raises(RuntimeError):
        msa.aligner_multiple([("a", "AC"), ("b", "AC")])


def test_moins_de_deux_sequences():
    with pytest.raises(ValueError):
        msa.aligner_multiple([("a", "AC")])


def test_programme_non_gere(tmp_path):
    with pytest.raises(ValueError):
        msa.aligner_multiple([("a", "AC"), ("b", "AC")], programme=str(tmp_path / "inconnu"))


def test_conservation():
    assert msa.conservation([("a", "AC-T"), ("b", "AGGT"), ("c", "ACGT")]) == "*  *"
    with pytest.raises(ValueError):
        msa.conservation([("a", "AC"), ("b", "A")])


def test_vers_fasta():
    assert msa.vers_fasta([("a", "AC-T"), ("b", "ACGT")], largeur=2) == ">a\nAC\n-T\n>b\nAC\nGT\n"
