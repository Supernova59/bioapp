from Bio.Align import PairwiseAligner


def aligner(seq_a, seq_b, mode="global"):
    """Aligne deux séquences. mode : 'global' ou 'local'.
    Retourne (score, ligne_a, ligne_b) avec des '-' pour les gaps."""
    if mode not in ("global", "local"):
        raise ValueError("mode doit être 'global' ou 'local'")
    al = PairwiseAligner()
    al.mode = mode
    al.match_score = 1
    al.mismatch_score = -1
    al.open_gap_score = -2
    al.extend_gap_score = -1
    meilleur = al.align(seq_a, seq_b)[0]

    blocs = list(zip(*meilleur.aligned))
    pa, pb = (blocs[0][0][0], blocs[0][1][0]) if mode == "local" else (0, 0)
    ligne_a, ligne_b = "", ""
    for (a0, a1), (b0, b1) in blocs:
        ligne_a += seq_a[pa:a0] + "-" * (b0 - pb) + seq_a[a0:a1]
        ligne_b += "-" * (a0 - pa) + seq_b[pb:b0] + seq_b[b0:b1]
        pa, pb = a1, b1
    if mode == "global":
        ligne_a += seq_a[pa:] + "-" * (len(seq_b) - pb)
        ligne_b += "-" * (len(seq_a) - pa) + seq_b[pb:]
    return meilleur.score, ligne_a, ligne_b
