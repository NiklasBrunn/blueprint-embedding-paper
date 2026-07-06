import subprocess, json
from pathlib import Path
REPO = Path(__file__).resolve().parents[2]
SH = REPO / "webapp" / "pipeline_driver.sh"

def _run(args):
    return subprocess.run(["bash", str(SH), *args, "--print-config"],
                          capture_output=True, text=True)

def test_scenario_input_workroot_override(tmp_path):
    scen = tmp_path / "scenario_config.jsonc"
    scen.write_text('{"species":["mouse"],"tissue":["cortex"],"celltype":null}')
    h5 = tmp_path / "data.h5ad"; h5.write_bytes(b"x")
    r = _run(["--model","gpt-5.5","--reasoning-effort","high",
              "--scenario-file",str(scen),"--input",str(h5),
              "--work-root",str(tmp_path/"runs"),
              "--non-interactive","--n-parallel","1"])
    assert r.returncode == 0, r.stderr
    cfg = json.loads(r.stdout)
    assert cfg["scenario_file"] == str(scen)
    assert cfg["blueprint_input"] == str(h5)
    assert cfg["work_root"] == str(tmp_path/"runs")
    assert cfg["n_parallel"] == 1
