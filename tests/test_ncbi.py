import io
import pytest
from Bio import Entrez
from core import ncbi


@pytest.fixture
def appels(monkeypatch):
    liste = []

    def faux_efetch(**kw):
        liste.append(kw)
        return io.StringIO(">NM_1 test\nACGT\n")

    monkeypatch.setattr(Entrez, "efetch", faux_efetch)
    return liste


def test_telechargement_puis_cache(tmp_path, appels):
    p1 = ncbi.telecharger_fasta("NM_1", "a@b.fr", cache=tmp_path)
    p2 = ncbi.telecharger_fasta("NM_1", "a@b.fr", cache=tmp_path)
    assert p1 == p2 and p1.read_text().startswith(">NM_1")
    assert len(appels) == 1  # le second appel vient du cache


def test_base_proteine(tmp_path, appels):
    ncbi.telecharger_fasta("NP_000537", "a@b.fr", cache=tmp_path)
    assert appels[0]["db"] == "protein"


def test_accession_invalide(tmp_path):
    with pytest.raises(ValueError):
        ncbi.telecharger_fasta("../etc/passwd", "a@b.fr", cache=tmp_path)


def test_email_requis(tmp_path, appels):
    with pytest.raises(ValueError):
        ncbi.telecharger_fasta("NM_1", "", cache=tmp_path)


class _Poignee:
    def __init__(self, donnees):
        self.donnees = donnees

    def close(self):
        pass


def _faux_entrez(monkeypatch, ids, resumes):
    appels = []
    monkeypatch.setattr(Entrez, "esearch", lambda **kw: appels.append(("esearch", kw)) or _Poignee({"IdList": ids}))
    monkeypatch.setattr(Entrez, "esummary", lambda **kw: appels.append(("esummary", kw)) or _Poignee(resumes))
    monkeypatch.setattr(Entrez, "read", lambda h: h.donnees)
    return appels


def test_rechercher(monkeypatch):
    appels = _faux_entrez(monkeypatch, ["11", "22"], [
        {"Id": "11", "AccessionVersion": "NM_1.1", "Title": "gène un", "Length": "100"},
        {"Id": "22", "Caption": "XM_2", "Title": "gène deux", "Length": "250"}])
    r = ncbi.rechercher("TP53[Gene]", "a@b.fr", base="nucleotide", maximum=5)
    assert r == [{"id": "11", "accession": "NM_1.1", "titre": "gène un", "longueur": 100},
                 {"id": "22", "accession": "XM_2", "titre": "gène deux", "longueur": 250}]
    assert appels[0][1]["retmax"] == 5 and appels[1][1]["id"] == "11,22"


def test_rechercher_sans_resultat(monkeypatch):
    appels = _faux_entrez(monkeypatch, [], [])
    assert ncbi.rechercher("introuvable", "a@b.fr") == []
    assert [a[0] for a in appels] == ["esearch"]


def test_rechercher_arguments_invalides():
    for args in (("", "a@b.fr"), ("x", ""), ("x", "a@b.fr", "genome")):
        with pytest.raises(ValueError):
            ncbi.rechercher(*args)


def test_base_imposee(tmp_path, appels):
    ncbi.telecharger_fasta("P10599", "a@b.fr", cache=tmp_path, base="protein")
    assert appels[0]["db"] == "protein"
    with pytest.raises(ValueError):
        ncbi.telecharger_fasta("P10599", "a@b.fr", cache=tmp_path / "x", base="genome")
