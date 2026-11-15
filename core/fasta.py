from Bio import SeqIO

def lire_fasta(chemin):
    """Retourne une liste de (identifiant, séquence) pour chaque entrée du fichier."""
    return [(r.id, str(r.seq)) for r in SeqIO.parse(chemin, "fasta")]