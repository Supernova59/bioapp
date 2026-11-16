from itertools import combinations

from core.alignement import aligner_complet, type_sequence


def pourcentage_identite(ligne_a, ligne_b):
    """Pourcentage de colonnes identiques (hors gaps) sur la longueur de l'alignement."""
    if not ligne_a:
        return 0.0
    ident = sum(1 for x, y in zip(ligne_a, ligne_b) if x == y and x != "-")
    return 100.0 * ident / len(ligne_a)


def choisir_bareme(seq_a, seq_b, match=1, mismatch=-1, ouverture=-2, extension=-1):
    """Retourne (paramètres pour aligner_complet, description).
    Protéines : BLOSUM62 avec gap -10/-1 ; sinon barème nucléotidique fourni."""
    if "proteine" in (type_sequence(seq_a), type_sequence(seq_b)):
        return ({"ouverture": -10, "extension": -1, "matrice": "BLOSUM62"},
                "protéines : BLOSUM62, gap -10/-1")
    return ({"match": match, "mismatch": mismatch,
             "ouverture": ouverture, "extension": extension},
            f"ADN : match {match}, mismatch {mismatch}, gap {ouverture}/{extension}")


def comparer(sequences, mode="global", **bareme):
    """Alignements deux à deux de 2 ou 3 séquences.
    sequences : dict {nom: séquence}. Retourne une liste de dicts, un par paire."""
    if not 2 <= len(sequences) <= 3:
        raise ValueError("il faut 2 ou 3 séquences")
    resultats = []
    for (nom_a, a), (nom_b, b) in combinations(sequences.items(), 2):
        a, b = a.upper(), b.upper()
        params, description = choisir_bareme(a, b, **bareme)
        score, la, lb, pos_a, pos_b = aligner_complet(a, b, mode, **params)
        resultats.append({
            "a": nom_a, "b": nom_b, "score": score,
            "ligne_a": la, "ligne_b": lb, "pos_a": pos_a, "pos_b": pos_b,
            "identite": pourcentage_identite(la, lb), "bareme": description,
        })
    return resultats
