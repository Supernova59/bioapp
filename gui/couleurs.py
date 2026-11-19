from html import escape

VERT, ROUGE, ORANGE, GRIS = "#2e9e4f", "#d6453d", "#d98a00", "#888888"


def milieu_alignement(ligne_a, ligne_b):
    """Ligne centrale : '|' pour une identité, espace sinon."""
    return "".join("|" if x == y and x != "-" else " " for x, y in zip(ligne_a, ligne_b))


def blocs_alignement(ligne_a, ligne_b, largeur=60):
    """Lignes de texte de l'alignement, par blocs : A, milieu, B, ligne vide."""
    milieu = milieu_alignement(ligne_a, ligne_b)
    lignes = []
    for i in range(0, len(ligne_a), largeur):
        lignes += [ligne_a[i:i + largeur], milieu[i:i + largeur], ligne_b[i:i + largeur], ""]
    return lignes


def _span(texte, couleur):
    return f'<span style="color:{couleur}">{escape(texte)}</span>'


def _ligne_sequence(sequence, milieu):
    morceaux = []
    for c, m in zip(sequence, milieu):
        morceaux.append(_span(c, ORANGE if c == "-" else VERT if m == "|" else ROUGE))
    return "".join(morceaux)


def colorer_alignement(lignes):
    """HTML coloré : identités en vert, différences en rouge, gaps en orange.
    `lignes` : en-tête, une ligne vide, puis des blocs (A, milieu, B, ligne vide)."""
    i = lignes.index("") if "" in lignes else len(lignes)
    sortie = [escape(ligne) for ligne in lignes[:i]]
    reste = lignes[i + 1:]
    for k in range(0, len(reste), 4):
        bloc = reste[k:k + 4]
        if len(bloc) < 3:
            sortie += [escape(ligne) for ligne in bloc]
            continue
        la, mid, lb = bloc[:3]
        sortie += ["", _ligne_sequence(la, mid), _span(mid, GRIS), _ligne_sequence(lb, mid)]
    return '<pre style="font-family:monospace; margin:0">' + "\n".join(sortie) + "</pre>"
