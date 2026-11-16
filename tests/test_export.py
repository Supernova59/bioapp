from core.export import fasta_aligne
from core.fasta import lire_fasta


def test_fasta_aligne_simple():
    assert fasta_aligne("a", "AC-T", "b", "ACGT") == ">a\nAC-T\n>b\nACGT\n"


def test_fasta_aligne_largeur():
    assert fasta_aligne("a", "ACGT", "b", "ACGT", largeur=2) == ">a\nAC\nGT\n>b\nAC\nGT\n"


def test_relecture_avec_gaps(tmp_path):
    f = tmp_path / "x.fasta"
    f.write_text(fasta_aligne("a", "AC-T", "b", "ACGT"))
    assert lire_fasta(f) == [("a", "AC-T"), ("b", "ACGT")]
