from pathlib import Path
from Bio import Entrez

CACHE = Path.home() / ".cache" / "bioapp"
PREFIXES_PROTEINE = ("NP_", "XP_", "YP_", "WP_")


def telecharger_fasta(accession, email, cache=CACHE):
    """Télécharge un FASTA depuis NCBI (avec cache disque) et retourne son chemin."""
    accession = accession.strip()
    if not accession or not accession.replace("_", "").replace(".", "").isalnum():
        raise ValueError("accession invalide")
    cache = Path(cache)
    cache.mkdir(parents=True, exist_ok=True)
    chemin = cache / f"{accession}.fasta"
    if chemin.exists():
        return chemin
    if not email:
        raise ValueError("NCBI demande une adresse e-mail")
    Entrez.email = email
    db = "protein" if accession.startswith(PREFIXES_PROTEINE) else "nucleotide"
    handle = Entrez.efetch(db=db, id=accession, rettype="fasta", retmode="text")
    try:
        texte = handle.read()
    finally:
        handle.close()
    if not texte.startswith(">"):
        raise ValueError(f"réponse NCBI inattendue pour {accession}")
    chemin.write_text(texte)
    return chemin
