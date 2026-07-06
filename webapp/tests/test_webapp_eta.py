import json
from webapp import eta


def test_known_combo_returns_estimate(tmp_path, monkeypatch):
    table = {"gpt_5_5|high|0": {"stage_a_s": 300.0, "stage_b_s": None, "n": 4}}
    p = tmp_path / "eta_table.json"
    p.write_text(json.dumps(table))
    monkeypatch.setattr(eta, "_TABLE_PATH", p)
    eta._load.cache_clear()
    seconds, n = eta.estimate("gpt-5.5", "high", refine=False)
    assert seconds == 300.0 and n == 4


def test_refine_adds_stage_b(tmp_path, monkeypatch):
    table = {"gpt_5_5|high|1": {"stage_a_s": 300.0, "stage_b_s": 120.0, "n": 3}}
    p = tmp_path / "eta_table.json"
    p.write_text(json.dumps(table))
    monkeypatch.setattr(eta, "_TABLE_PATH", p)
    eta._load.cache_clear()
    seconds, n = eta.estimate("gpt-5.5", "high", refine=True)
    assert seconds == 420.0 and n == 3


def test_missing_combo_returns_none(tmp_path, monkeypatch):
    p = tmp_path / "eta_table.json"
    p.write_text("{}")
    monkeypatch.setattr(eta, "_TABLE_PATH", p)
    eta._load.cache_clear()
    seconds, n = eta.estimate("gpt-9", "low", refine=False)
    assert seconds is None and n == 0


def test_format_eta_human():
    assert "min" in eta.format_eta(300.0, 4)
    assert eta.format_eta(None, 0) == "No timing data for this combination yet."
