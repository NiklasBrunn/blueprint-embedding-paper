"""Native 'Browse…' file chooser for the local-first app.

Browsers never expose a selected file's absolute path (security), so a real
"pick in Finder → read from disk" requires the *server* to open a native dialog.
That only makes sense because the app is local-first: the server runs on the
same machine as the browser. On macOS we shell out to ``osascript`` (no GUI
dependency, runs out-of-process so it is safe from a request thread); elsewhere
we fall back to a tkinter dialog if available.
"""
from __future__ import annotations

import platform
import subprocess


def pick_file_native(prompt: str = "Select an .h5ad file") -> str:
    """Open a native file chooser on the server's machine and return the chosen
    absolute path, or '' if cancelled or unavailable. Blocks until the user
    picks (deliberate user action)."""
    if platform.system() == "Darwin":
        return _pick_macos(prompt)
    return _pick_tk(prompt)


def _pick_macos(prompt: str) -> str:
    script = f'POSIX path of (choose file with prompt "{prompt}")'
    try:
        out = subprocess.run(["osascript", "-e", script],
                             capture_output=True, text=True, timeout=600)
    except Exception:
        return ""
    # rc != 0 on cancel ("User canceled. (-128)") -> no selection.
    return out.stdout.strip() if out.returncode == 0 else ""


def _pick_tk(prompt: str) -> str:
    try:
        import tkinter as tk
        from tkinter import filedialog
    except Exception:
        return ""
    try:
        root = tk.Tk()
        root.withdraw()
        path = filedialog.askopenfilename(
            title=prompt,
            filetypes=[("AnnData h5ad", "*.h5ad *.h5ad.gz"), ("All files", "*.*")])
        root.destroy()
        return path or ""
    except Exception:
        return ""
