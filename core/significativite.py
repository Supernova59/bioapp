import random
from statistics import mean, pstdev

from Bio.Align import PairwiseAligner, substitution_matrices


def _aligneur(mode, match, mismatch, ouverture, extension, matrice):
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
    return al


def significativite(seq_a, seq_b, melanges=200, graine=None, mode="local",
                    match=1, mismatch=-1, ouverture=-2, extension=-1, matrice=None):
    """Test par mélange (comme PRSS) : la séquence B est mélangée `melanges` fois et
    réalignée avec A ; on compare le score réel à ceux des séquences mélangées.
    Retourne un dict : score, moyenne, ecart_type, z, meilleurs, melanges, p (empirique)."""
    if melanges < 1:
        raise ValueError("melanges doit être >= 1")
    seq_a, seq_b = seq_a.upper(), seq_b.upper()
    al = _aligneur(mode, match, mismatch, ouverture, extension, matrice)
    reel = al.score(seq_a, seq_b)
    rng = random.Random(graine)
    lettres = list(seq_b)
    scores = []
    for _ in range(melanges):
        rng.shuffle(lettres)
        scores.append(al.score(seq_a, "".join(lettres)))
    moy, sd = mean(scores), pstdev(scores)
    meilleurs = sum(1 for s in scores if s >= reel)
    return {
        "score": reel, "moyenne": moy, "ecart_type": sd,
        "z": (reel - moy) / sd if sd > 0 else (float("inf") if reel > moy else 0.0),
        "meilleurs": meilleurs, "melanges": melanges,
        "p": (meilleurs + 1) / (melanges + 1),
    }
