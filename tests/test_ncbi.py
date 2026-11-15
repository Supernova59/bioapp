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
