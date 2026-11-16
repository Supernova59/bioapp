import re

from core.sequences import complement_inverse

_ELEMENT = re.compile(r"(x|[A-Z]|\[[A-Z>]+\]|\{[A-Z]+\})(?:\((\d+)(?:,(\d+))?\))?")


def _element(el):
    prefixe = suffixe = ""
    if el.startswith("<"):
        prefixe, el = "^", el[1:]
    if el.endswith(">"):
        suffixe, el = "$", el[:-1]
    if el == "":
        return prefixe + suffixe
    m = _ELEMENT.fullmatch(el)
    if not m:
        raise ValueError(f"élément Prosite invalide : {el!r}")
    base, n, p = m.groups()
    if base == "x":
        rx = "."
    elif base[0] == "[":
        lettres = base[1:-1]
        if ">" in lettres:
            seuls = lettres.replace(">", "")
            rx = f"(?:[{seuls}]|$)" if seuls else "$"
        else:
            rx = base
    elif base[0] == "{":
        rx = "[^" + base[1:-1] + "]"
    else:
        rx = base
    rep = "" if n is None else (f"{{{n}}}" if p is None else f"{{{n},{p}}}")
    return prefixe + rx + rep + suffixe


def prosite_vers_regex(motif):
    """Traduit un motif au format Prosite en expression régulière Python.
    Ex. 'C-x(2,4)-[LIV]-{P}-G.' -> 'C.{2,4}[LIV][^P]G'"""
    motif = motif.strip().rstrip(".")
    if not motif:
        raise ValueError("motif vide")
    return "".join(_element(el) for el in motif.split("-"))


def _recherche(rx, seq, brin):
    trouves = []
    for m in rx.finditer(seq):
        sous = m.group(1)
        trouves.append({"debut": m.start() + 1, "fin": m.start() + len(sous),
                        "sequence": sous, "brin": brin})
    return trouves


def chercher_motif(sequence, motif, deux_brins=False):
    """Cherche toutes les occurrences (y compris chevauchantes) d'un motif Prosite.
    Positions en base 1, fin incluse. deux_brins : cherche aussi le complément inverse
    (ADN) ; les positions du brin '-' sont exprimées sur la séquence d'origine."""
    rx = re.compile("(?=(" + prosite_vers_regex(motif) + "))")
    seq = sequence.upper()
    resultats = _recherche(rx, seq, "+")
    if deux_brins:
        n = len(seq)
        for r in _recherche(rx, complement_inverse(seq), "-"):
            r["debut"], r["fin"] = n - r["fin"] + 1, n - r["debut"] + 1
            resultats.append(r)
    return sorted(resultats, key=lambda r: (r["debut"], r["brin"]))
