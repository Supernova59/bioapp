import os
import shutil
import subprocess
import tempfile

from core.fasta import lire_fasta

PROGRAMMES = ("clustalo", "mafft")


def trouver_programme():
    """Chemin du premier programme d'alignement multiple installé, ou None."""
    for nom in PROGRAMMES:
        chemin = shutil.which(nom)
        if chemin:
            return chemin
    return None


def _commande(programme, entree, sortie):
    nom = os.path.basename(programme)
    if nom.startswith("clustalo"):
        return [programme, "-i", entree, "-o", sortie, "--outfmt=fasta", "--force"], None
    if nom.startswith("mafft"):
        return [programme, "--auto", "--quiet", entree], sortie   # MAFFT écrit sur stdout
    raise ValueError(f"programme non géré : {nom}")


def aligner_multiple(sequences, programme=None, delai=300):
    """Alignement multiple via Clustal Omega ou MAFFT (installés sur la machine).
    sequences : liste de (identifiant, séquence), au moins 2.
    Retourne une liste de (identifiant, séquence alignée avec '-')."""
    if len(sequences) < 2:
        raise ValueError("il faut au moins 2 séquences")
    programme = programme or trouver_programme()
    if programme is None:
        raise RuntimeError("aucun programme d'alignement multiple trouvé "
                           "(installe-le avec : sudo apt install clustalo)")
    with tempfile.TemporaryDirectory() as dossier:
        entree = os.path.join(dossier, "entree.fasta")
        sortie = os.path.join(dossier, "sortie.fasta")
        with open(entree, "w") as f:
            for i, (ident, seq) in enumerate(sequences):
                f.write(f">s{i}\n{seq}\n")   # identifiants neutres : pas de caractère spécial
        commande, fichier_stdout = _commande(programme, entree, sortie)
        res = subprocess.run(commande, capture_output=True, text=True, timeout=delai)
        if res.returncode != 0:
            raise RuntimeError(f"{os.path.basename(programme)} a échoué : {res.stderr.strip()[:300]}")
        if fichier_stdout:
            with open(fichier_stdout, "w") as f:
                f.write(res.stdout)
        resultat = {ident: seq for ident, seq in lire_fasta(sortie)}
    try:
        return [(ident, resultat[f"s{i}"]) for i, (ident, _) in enumerate(sequences)]
    except KeyError:
        raise RuntimeError("sortie inattendue du programme d'alignement")


def conservation(alignement):
    """Ligne de conservation : '*' si toutes les séquences sont identiques et sans gap."""
    lignes = [seq for _, seq in alignement]
    if len({len(s) for s in lignes}) != 1:
        raise ValueError("les séquences alignées n'ont pas la même longueur")
    return "".join("*" if len(set(col)) == 1 and "-" not in col else " "
                   for col in zip(*lignes))


def vers_fasta(alignement, largeur=60):
    """Texte FASTA d'un alignement multiple (les '-' sont conservés)."""
    sortie = []
    for ident, seq in alignement:
        sortie.append(f">{ident}")
        sortie += [seq[i:i + largeur] for i in range(0, len(seq), largeur)]
    return "\n".join(sortie) + "\n"
