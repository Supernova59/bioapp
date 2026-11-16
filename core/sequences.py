from Bio.Seq import Seq


def complement_inverse(seq):
    """Complément inverse d'une séquence nucléotidique (la casse est conservée)."""
    return str(Seq(seq).reverse_complement())
