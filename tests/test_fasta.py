from core.fasta import lire_fasta

def test_lire_fasta(tmp_path):
    f = tmp_path / "x.fasta"
    f.write_text(">seq1\nACGT\nACGT\n>seq2\nTTTT\n")
    assert lire_fasta(f) == [("seq1", "ACGTACGT"), ("seq2", "TTTT")]