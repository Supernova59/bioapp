import io

from Bio.Blast import NCBIWWW, NCBIXML

from core.alignement import type_sequence

# programme -> bases proposées (la première est la valeur par défaut)
PROGRAMMES = {
    "blastn": ("nt", "refseq_rna"),
    "blastp": ("swissprot", "nr", "refseq_protein", "pdb"),
    "blastx": ("swissprot", "nr", "refseq_protein"),
}
LONGUEUR_MAX = 50_000


def analyser_xml(source, maximum=10):
    """Lit un résultat BLAST au format XML (fichier ou flux texte).
    Retourne une liste de hits (le meilleur HSP de chacun), du meilleur au moins bon."""
    try:
        record = NCBIXML.read(source)
    except Exception as e:
        raise RuntimeError(f"résultat BLAST illisible : {e}")
    longueur = record.query_length or 1
    hits = []
    for aln in record.alignments[:maximum]:
        hsp = aln.hsps[0]
        hits.append({
            "accession": aln.accession,
            "titre": aln.hit_def or aln.title,
            "longueur": aln.length,
            "evalue": hsp.expect,
            "score": hsp.score,
            "bits": hsp.bits,
            "identite": 100.0 * hsp.identities / hsp.align_length,
            "couverture": 100.0 * (abs(hsp.query_end - hsp.query_start) + 1) / longueur,
            "debut_requete": hsp.query_start, "fin_requete": hsp.query_end,
            "debut_cible": hsp.sbjct_start, "fin_cible": hsp.sbjct_end,
            "requete": hsp.query, "milieu": hsp.match, "cible": hsp.sbjct,
        })
    return hits


def lancer_blast(sequence, programme="blastp", base=None, e_value=0.01, maximum=10):
    """Envoie la séquence au BLAST web de NCBI (peut prendre plusieurs minutes)."""
    if programme not in PROGRAMMES:
        raise ValueError("programme doit être blastn, blastp ou blastx")
    base = base or PROGRAMMES[programme][0]
    if base not in PROGRAMMES[programme]:
        raise ValueError(f"base {base!r} non proposée pour {programme}")
    sequence = "".join(sequence.split()).upper()
    if len(sequence) < 10:
        raise ValueError("séquence trop courte (au moins 10)")
    if len(sequence) > LONGUEUR_MAX:
        raise ValueError(f"séquence trop longue (maximum {LONGUEUR_MAX})")
    adn = type_sequence(sequence) == "nucleotide"
    if programme in ("blastn", "blastx") and not adn:
        raise ValueError(f"{programme} demande une séquence d'ADN")
    if programme == "blastp" and adn:
        raise ValueError("blastp demande une protéine")
    handle = NCBIWWW.qblast(programme, base, sequence, expect=e_value,
                            hitlist_size=maximum, format_type="XML")
    try:
        texte = handle.read()
    finally:
        handle.close()
    return analyser_xml(io.StringIO(texte), maximum)
