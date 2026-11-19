from Bio.Seq import Seq
from Bio.SeqUtils.ProtParamData import kd

from core.sequences import complement_inverse


def profil_hydropathie(proteine, fenetre=9):
    """Profil de Kyte-Doolittle : liste de (position centrale en base 1, score moyen)."""
    proteine = proteine.upper().replace("*", "")
    if fenetre < 1 or fenetre % 2 == 0:
        raise ValueError("la fenêtre doit être un entier impair")
    if len(proteine) < fenetre:
        raise ValueError("fenêtre plus grande que la séquence")
    inconnus = sorted(set(proteine) - set(kd))
    if inconnus:
        raise ValueError(f"acide(s) aminé(s) inconnu(s) : {''.join(inconnus)}")
    valeurs = [kd[a] for a in proteine]
    return [(i + fenetre // 2 + 1, sum(valeurs[i:i + fenetre]) / fenetre)
            for i in range(len(valeurs) - fenetre + 1)]


def regions_hydrophobes(proteine, fenetre=19, seuil=1.6):
    """Régions candidates (hélices transmembranaires : fenêtre 19, seuil 1,6).
    Retourne une liste de (debut, fin) en base 1, fin incluse."""
    profil = profil_hydropathie(proteine, fenetre)
    demi = fenetre // 2
    regions, debut, fin = [], None, None
    for centre, score in profil:
        if score >= seuil:
            if debut is None:
                debut = centre - demi
            fin = centre + demi
        elif debut is not None:
            regions.append((debut, fin))
            debut = None
    if debut is not None:
        regions.append((debut, fin))
    return regions


def traduire(adn, cadre=1, table=1):
    """Traduit un cadre de lecture (1, 2, 3 sur le brin +, -1, -2, -3 sur le brin -).
    Le dernier codon incomplet est ignoré ; '*' représente un codon stop."""
    if cadre not in (1, 2, 3, -1, -2, -3):
        raise ValueError("cadre doit valoir 1, 2, 3, -1, -2 ou -3")
    seq = adn.upper().replace("U", "T")
    if cadre < 0:
        seq = complement_inverse(seq)
    seq = seq[abs(cadre) - 1:]
    seq = seq[:len(seq) - len(seq) % 3]
    return str(Seq(seq).translate(table=table))


def chercher_orfs(adn, longueur_min=30, table=1, deux_brins=True):
    """ORF complets (de ATG jusqu'à un stop inclus) dans les 3 ou 6 cadres.
    Positions en base 1 sur la séquence d'origine, fin incluse (stop compris).
    Retourne une liste de dicts triée par position."""
    adn = adn.upper().replace("U", "T")
    n = len(adn)
    cadres = (1, 2, 3, -1, -2, -3) if deux_brins else (1, 2, 3)
    orfs = []
    for cadre in cadres:
        prot = traduire(adn, cadre, table)
        decalage = abs(cadre) - 1
        i = 0
        while i < len(prot):
            j = prot.find("*", i)
            if j == -1:
                break               # pas de stop : ORF incomplet, ignoré
            m = prot.find("M", i, j)
            if m != -1 and j - m >= longueur_min:
                debut = decalage + 3 * m + 1
                fin = decalage + 3 * (j + 1)
                if cadre < 0:
                    debut, fin = n - fin + 1, n - debut + 1
                orfs.append({"brin": "+" if cadre > 0 else "-", "cadre": cadre,
                             "debut": debut, "fin": fin, "longueur_aa": j - m,
                             "proteine": prot[m:j]})
            i = j + 1
    return sorted(orfs, key=lambda o: (o["debut"], o["brin"]))
