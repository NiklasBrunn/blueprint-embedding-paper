from webapp.jobs import JobRegistry

def test_job_runs_and_completes(tmp_path):
    reg = JobRegistry(); log = tmp_path/"job.log"
    jid = reg.start_subprocess(["bash","-c","echo hello"], cwd=tmp_path,
                               env=None, log_path=log)
    reg.wait(jid, timeout=10)
    assert reg.get(jid).status == "done" and "hello" in log.read_text()

def test_tail_recent_lines(tmp_path):
    reg = JobRegistry(); log = tmp_path/"job.log"
    jid = reg.start_subprocess(["bash","-c","printf 'a\\nb\\nc\\n'"],
                               cwd=tmp_path, env=None, log_path=log)
    reg.wait(jid, timeout=10)
    assert reg.tail(jid, n=2) == ["b", "c"]

def test_failed_job_marked_failed(tmp_path):
    reg = JobRegistry(); log = tmp_path/"job.log"
    jid = reg.start_subprocess(["bash","-c","exit 3"], cwd=tmp_path,
                               env=None, log_path=log)
    reg.wait(jid, timeout=10)
    assert reg.get(jid).status == "failed"
