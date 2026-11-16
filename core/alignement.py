from Bio.Align import PairwiseAligner, substitution_matrices

VIDE = (0.0, "", "", (0, 0), (0, 0))
MODES = ("global", "local", "semi-global")


def type_sequence(seq):
    """'nucleotide' si la séquence ne contient que A, C, G, T, U, N, sinon 'proteine'."""
    lettres = set(seq.upper()) - set("-*")
    return "nucleotide" if lettres <= set("ACGTUN") else "proteine"


def aligner_complet(seq_a, seq_b, mode="global", match=1, mismatch=-1,
                    ouverture=-2, extension=-1, matrice=None):
    """Aligne deux séquences. mode : 'global', 'local' ou 'semi-global'
    (semi-global : les gaps aux extrémités ne sont pas pénalisés).
    matrice : None (barème match/mismatch) ou nom d'une matrice, par ex. 'BLOSUM62'.
    Retourne (score, ligne_a, ligne_b, (debut_a, fin_a), (debut_b, fin_b)).
    Positions en base 0, fin exclue, sur la région alignée."""
    if mode not in MODES:
        raise ValueError("mode doit être 'global', 'local' ou 'semi-global'")
    al = PairwiseAligner()
    al.mode = "local" if mode == "local" else "global"
    if matrice:
        al.substitution_matrix = substitution_matrices.load(matrice)
    else:
        al.match_score = match
        al.mismatch_score = mismatch
    al.open_gap_score = ouverture
    al.extend_gap_score = extension
    if mode == "semi-global":
        al.end_gap_score = 0.0
    try:
        meilleur = al.align(seq_a, seq_b)[0]
    except (IndexError, ValueError):
        return VIDE
    blocs = list(zip(*meilleur.aligned))
    if not blocs:
        return VIDE

    deb_a, deb_b = int(blocs[0][0][0]), int(blocs[0][1][0])
    fin_a, fin_b = int(blocs[-1][0][1]), int(blocs[-1][1][1])
    pa, pb = (deb_a, deb_b) if mode == "local" else (0, 0)
    ligne_a, ligne_b = "", ""
    for (a0, a1), (b0, b1) in blocs:
        ligne_a += seq_a[pa:a0] + "-" * (b0 - pb) + seq_a[a0:a1]
        ligne_b += "-" * (a0 - pa) + seq_b[pb:b0] + seq_b[b0:b1]
        pa, pb = a1, b1
    if mode != "local":
        ligne_a += seq_a[pa:] + "-" * (len(seq_b) - pb)
        ligne_b += "-" * (len(seq_a) - pa) + seq_b[pb:]
        return meilleur.score, ligne_a, ligne_b, (0, len(seq_a)), (0, len(seq_b))
    return meilleur.score, ligne_a, ligne_b, (deb_a, fin_a), (deb_b, fin_b)


def aligner(seq_a, seq_b, mode="global"):
    """Version simple : (score, ligne_a, ligne_b) avec le barème par défaut."""
    return aligner_complet(seq_a, seq_b, mode)[:3]
