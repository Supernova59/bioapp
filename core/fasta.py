from Bio import SeqIO


def lire_fasta(chemin):
    """Retourne une liste de (identifiant, séquence) pour chaque entrée du fichier."""
    with open(chemin) as f:
        return [(r.id, str(r.seq)) for r in SeqIO.parse(f, "fasta")]
