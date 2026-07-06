"""In-memory background job registry for long-running subprocesses."""
from __future__ import annotations
import subprocess, threading, time, uuid
from dataclasses import dataclass, field
from pathlib import Path

@dataclass
class Job:
    id: str
    log_path: Path
    proc: "subprocess.Popen | None" = None
    status: str = "running"   # running | done | failed | cancelled
    meta: dict = field(default_factory=dict)

class JobRegistry:
    def __init__(self):
        self._jobs: dict[str, Job] = {}
        self._lock = threading.Lock()

    def start_subprocess(self, cmd, *, cwd, env, log_path, meta=None) -> str:
        jid = uuid.uuid4().hex[:12]
        job = Job(id=jid, log_path=Path(log_path), meta=meta or {})
        with self._lock:
            self._jobs[jid] = job
        def _run():
            with open(log_path, "w") as fh:
                job.proc = subprocess.Popen(cmd, cwd=str(cwd), env=env,
                                            stdout=fh, stderr=subprocess.STDOUT,
                                            text=True)
                rc = job.proc.wait()
            # A cancel() may have already set status; don't clobber it.
            if job.status == "running":
                job.status = "done" if rc == 0 else "failed"
        threading.Thread(target=_run, daemon=True).start()
        return jid

    def get(self, jid) -> Job:
        return self._jobs[jid]

    def tail(self, jid, n=40) -> list:
        p = self._jobs[jid].log_path
        if not p.exists():
            return []
        return p.read_text().splitlines()[-n:]

    def cancel(self, jid):
        job = self._jobs[jid]
        if job.proc is not None and job.status == "running":
            job.proc.terminate()
            job.status = "cancelled"

    def wait(self, jid, timeout=None):
        t0 = time.time()
        while self._jobs[jid].status == "running":
            if timeout and time.time() - t0 > timeout:
                raise TimeoutError(jid)
            time.sleep(0.05)
