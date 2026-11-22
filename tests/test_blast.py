import io
import pytest
from core import blast


def _hsp(qf, qt, hf, ht, ident, alen, ev, bits=100.0, score=200, qseq="MSHY", hseq="MSHY", mid="MSHY"):
    return (f"<Hsp><Hsp_num>1</Hsp_num><Hsp_bit-score>{bits}</Hsp_bit-score><Hsp_score>{score}</Hsp_score>"
            f"<Hsp_evalue>{ev}</Hsp_evalue><Hsp_query-from>{qf}</Hsp_query-from><Hsp_query-to>{qt}</Hsp_query-to>"
            f"<Hsp_hit-from>{hf}</Hsp_hit-from><Hsp_hit-to>{ht}</Hsp_hit-to><Hsp_query-frame>0</Hsp_query-frame>"
            f"<Hsp_hit-frame>0</Hsp_hit-frame><Hsp_identity>{ident}</Hsp_identity>"
            f"<Hsp_positive>{ident}</Hsp_positive><Hsp_gaps>0</Hsp_gaps><Hsp_align-len>{alen}</Hsp_align-len>"
            f"<Hsp_qseq>{qseq}</Hsp_qseq><Hsp_hseq>{hseq}</Hsp_hseq><Hsp_midline>{mid}</Hsp_midline></Hsp>")


def _hit(num, acc, titre, longueur, hsp):
    return (f"<Hit><Hit_num>{num}</Hit_num><Hit_id>{acc}</Hit_id><Hit_def>{titre}</Hit_def>"
            f"<Hit_accession>{acc}</Hit_accession><Hit_len>{longueur}</Hit_len><Hit_hsps>{hsp}</Hit_hsps></Hit>")


def _xml(hits=""):
    return ('<?xml version="1.0"?><BlastOutput><BlastOutput_program>blastp</BlastOutput_program>'
            '<BlastOutput_query-ID>Query_1</BlastOutput_query-ID><BlastOutput_query-def>q</BlastOutput_query-def>'
            '<BlastOutput_query-len>100</BlastOutput_query-len><BlastOutput_param><Parameters>'
            '<Parameters_expect>0.01</Parameters_expect></Parameters></BlastOutput_param>'
            '<BlastOutput_iterations><Iteration><Iteration_iter-num>1</Iteration_iter-num>'
            f'<Iteration_hits>{hits}</Iteration_hits></Iteration></BlastOutput_iterations></BlastOutput>')


DEUX_HITS = (_hit(1, "P10599", "Thioredoxin", 105, _hsp(1, 100, 5, 104, 80, 100, "2.5e-60"))
             + _hit(2, "P0AA25", "Thioredoxin 1", 109, _hsp(20, 69, 10, 59, 25, 50, "1e-10", mid="MS Y")))


def test_analyser_deux_hits():
    h = blast.analyser_xml(io.StringIO(_xml(DEUX_HITS)))
    assert [x["accession"] for x in h] == ["P10599", "P0AA25"]
    assert h[0]["evalue"] == 2.5e-60 and h[0]["identite"] == 80.0 and h[0]["couverture"] == 100.0
    assert h[1]["couverture"] == 50.0 and h[1]["debut_cible"] == 10 and h[1]["fin_cible"] == 59
    assert h[1]["titre"] == "Thioredoxin 1" and h[1]["longueur"] == 109


def test_analyser_maximum_et_vide():
    assert len(blast.analyser_xml(io.StringIO(_xml(DEUX_HITS)), maximum=1)) == 1
    assert blast.analyser_xml(io.StringIO(_xml())) == []


def test_xml_illisible():
    with pytest.raises(RuntimeError):
        blast.analyser_xml(io.StringIO("<html>erreur</html>"))


def test_lancer_blast_avec_faux_qblast(monkeypatch):
    appels = []

    def faux_qblast(*args, **kw):
        appels.append((args, kw))
        return io.StringIO(_xml(DEUX_HITS))

    monkeypatch.setattr(blast.NCBIWWW, "qblast", faux_qblast)
    h = blast.lancer_blast("MSHY IELTEE\nNFES", programme="blastp", e_value=0.001, maximum=5)
    assert len(h) == 2
    args, kw = appels[0]
    assert args == ("blastp", "swissprot", "MSHYIELTEENFES")
    assert kw["expect"] == 0.001 and kw["hitlist_size"] == 5 and kw["format_type"] == "XML"


def test_lancer_blast_validations():
    cas = [("MSHYIELTEENFES", "blastz", None), ("MSHYIELTEENFES", "blastp", "nt"),
           ("MSHY", "blastp", None), ("ACGTACGTACGT", "blastp", None),
           ("MSHYIELTEENFES", "blastn", None), ("MSHYIELTEENFES", "blastx", None)]
    for sequence, programme, base in cas:
        with pytest.raises(ValueError):
            blast.lancer_blast(sequence, programme, base)
