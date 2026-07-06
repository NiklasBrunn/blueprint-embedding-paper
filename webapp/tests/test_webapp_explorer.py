from webapp.explorer import branches_from_summary


def test_branches_from_summary_returns_list():
    summary = {"blueprint": {"axes": [
        {"name": "epi_root", "kind": "plain", "gate": None, "positive": "epiblast"},
        {"name": "epi_mid", "kind": "gated", "gate": "epiblast", "positive": "mid"},
    ]}}
    branches = branches_from_summary(summary)
    assert isinstance(branches, list)
    # epi_root.positive == epi_mid.gate == "epiblast" -> a 2-node branch
    assert branches and branches[0]["axes"] == ["epi_root", "epi_mid"]
