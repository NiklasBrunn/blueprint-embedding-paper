"""The native dialog itself can't run in CI; mock the osascript subprocess."""
import webapp.filepicker as F


class _Res:
    def __init__(self, returncode, stdout="", stderr=""):
        self.returncode = returncode
        self.stdout = stdout
        self.stderr = stderr


def test_macos_returns_chosen_path(monkeypatch):
    monkeypatch.setattr(F.platform, "system", lambda: "Darwin")
    monkeypatch.setattr(F.subprocess, "run",
                        lambda *a, **k: _Res(0, "/tmp/data.h5ad\n"))
    assert F.pick_file_native() == "/tmp/data.h5ad"


def test_macos_cancel_returns_empty(monkeypatch):
    monkeypatch.setattr(F.platform, "system", lambda: "Darwin")
    monkeypatch.setattr(F.subprocess, "run",
                        lambda *a, **k: _Res(1, "", "User canceled. (-128)"))
    assert F.pick_file_native() == ""


def test_macos_osascript_unavailable_returns_empty(monkeypatch):
    monkeypatch.setattr(F.platform, "system", lambda: "Darwin")

    def boom(*a, **k):
        raise FileNotFoundError("osascript")

    monkeypatch.setattr(F.subprocess, "run", boom)
    assert F.pick_file_native() == ""
