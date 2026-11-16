def fasta_aligne(id_a, ligne_a, id_b, ligne_b, largeur=60):
    """FASTA de deux séquences alignées (les '-' sont conservés)."""
    def bloc(ident, ligne):
        morceaux = [ligne[i:i + largeur] for i in range(0, len(ligne), largeur)]
        return f">{ident}\n" + "\n".join(morceaux) + "\n"
    return bloc(id_a, ligne_a) + bloc(id_b, ligne_b)
