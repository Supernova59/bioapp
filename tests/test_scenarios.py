from core.comparaison import comparer
from gui.scenarios import SCENARIOS


def test_scores_attendus_des_scenarios():
    for sc in SCENARIOS:
        if sc["action"] != "aligner":
            continue
        r = comparer({"A": sc["A"][1], "B": sc["B"][1]}, mode=sc["mode"], **sc["params"])[0]
        assert r["score"] == sc["score"], sc["titre"]


def test_scenarios_bien_formes():
    for sc in SCENARIOS:
        assert sc["mode"] in ("global", "semi-global", "local")
        assert sc["action"] in ("aligner", "dotplot")
        assert set(sc["params"]) == {"match", "mismatch", "ouverture", "extension"}
        assert sc["attendu"]
