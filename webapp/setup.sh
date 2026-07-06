#!/usr/bin/env bash
# Create webapp/.venv and install the app's complete dependency set.
# Idempotent: re-running just re-installs into the existing venv.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VENV="$HERE/.venv"

# ---------------------------------------------------------------------------
# Preferred path: uv.
# uv builds/manages its own *standalone* (relocatable) Python interpreters. The
# stdlib `venv` module produces a BROKEN venv from such an interpreter — its
# stdlib path stays pinned to the build-time "/install" prefix, so the venv's
# python crashes on startup with "ModuleNotFoundError: No module named
# 'encodings'". uv creates the venv correctly and installs without needing
# ensurepip/pip in the venv at all. If uv is on PATH, just use it.
# ---------------------------------------------------------------------------
if command -v uv >/dev/null 2>&1; then
  echo "Using uv ($(uv --version)) to create the venv and install dependencies"
  # Honor an explicit PYTHON= request (version like 3.12, or a path); otherwise
  # let uv pick/download a suitable interpreter. The app needs >=3.10.
  if [ -n "${PYTHON:-}" ]; then
    uv venv --python "$PYTHON" "$VENV"
  else
    uv venv --python ">=3.10" "$VENV"
  fi
  uv pip install --python "$VENV/bin/python" -r "$HERE/requirements.txt"
  echo
  echo "Done. Launch with:"
  echo "  $VENV/bin/python -m webapp"
  exit 0
fi

# ---------------------------------------------------------------------------
# Fallback path: stdlib venv (no uv available).
# ---------------------------------------------------------------------------
# Pick a modern Python (>=3.10). The app needs anndata>=0.11 to read .h5ad files
# written by recent scanpy/anndata (the "null" encoding for None) — that anndata
# is unavailable on Python 3.9. Override by setting PYTHON=... yourself.
if [ -z "${PYTHON:-}" ]; then
  for cand in python3.12 python3.11 python3.13 python3.10 python3; do
    if command -v "$cand" >/dev/null 2>&1; then PYTHON="$cand"; break; fi
  done
fi
PYTHON="${PYTHON:-python3}"
echo "Using interpreter: ${PYTHON} ($(${PYTHON} --version 2>&1))"

# Create the venv. Some Python builds ship a broken/absent ensurepip (e.g. a
# Debian/Ubuntu python without the matching python3.x-venv package), which makes
# a plain `python -m venv` abort. In that case create the venv without pip and
# bootstrap pip ourselves.
if [ ! -d "$VENV/bin" ]; then
  echo "Creating venv at $VENV"
  rm -rf "$VENV"
  if ! "$PYTHON" -m venv "$VENV" 2>/dev/null; then
    echo "  (ensurepip unavailable — creating without pip and bootstrapping it)"
    rm -rf "$VENV"
    "$PYTHON" -m venv --without-pip "$VENV"
  fi
fi

# Sanity-check that the venv's interpreter actually runs. A relocatable/standalone
# Python (uv- or pyenv-managed) yields a venv whose python cannot find its stdlib
# ("No module named 'encodings'"); bootstrapping pip into it is hopeless, so fail
# early with an actionable message instead of a cryptic crash.
if ! "$VENV/bin/python" -c "import encodings" >/dev/null 2>&1; then
  rm -rf "$VENV"
  cat >&2 <<EOF

ERROR: the virtual environment's Python does not work — its standard library
cannot be found. This happens when '${PYTHON}' is a relocatable/standalone build
(e.g. a uv- or pyenv-managed interpreter) that the stdlib 'venv' module cannot
build a working venv from.

Fix (pick one):
  • Install uv, then re-run this script (it will use uv automatically):
        curl -LsSf https://astral.sh/uv/install.sh | sh
        bash webapp/setup.sh
  • Use a normal system Python:
        sudo apt install python3.12-venv          # Debian/Ubuntu, if needed
        PYTHON=/usr/bin/python3 bash webapp/setup.sh
  • If you use conda, point at its interpreter:
        PYTHON="\$(command -v python)" bash webapp/setup.sh
EOF
  exit 1
fi

# Make sure the venv has pip (no-op if venv already created it).
if [ ! -x "$VENV/bin/pip" ] && [ ! -f "$VENV/bin/pip3" ]; then
  echo "Bootstrapping pip into the venv"
  if ! "$VENV/bin/python" -m ensurepip --upgrade --default-pip >/dev/null 2>&1; then
    echo "  (ensurepip failed — fetching get-pip.py)"
    GETPIP="$(mktemp)"
    curl -fsSL https://bootstrap.pypa.io/get-pip.py -o "$GETPIP" \
      || { echo "ERROR: could not download get-pip.py (no network?)."; exit 1; }
    "$VENV/bin/python" "$GETPIP"
    rm -f "$GETPIP"
  fi
fi

"$VENV/bin/python" -m pip install --upgrade pip
"$VENV/bin/python" -m pip install -r "$HERE/requirements.txt"
echo
echo "Done. Launch with:"
echo "  $VENV/bin/python -m webapp"
