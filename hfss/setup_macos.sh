#!/bin/bash
# Build openEMS from source and set up the repo's virtual environment on macOS.
# There is no macOS release of openEMS, so this compiles it (roughly 10-30 minutes).
# Needs Homebrew (https://brew.sh) and Xcode Command Line Tools. Safe to re-run.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PREFIX="${OPENEMS_INSTALL_PATH:-$HOME/opt/openEMS}"
SRC="$HOME/opt/src/openEMS-Project"
PYTHON="${PYTHON:-python3.13}"

command -v brew >/dev/null || { echo "Homebrew not found: install it from https://brew.sh first." >&2; exit 1; }
brew install cmake boost hdf5

"$PYTHON" -m venv "$HERE/../.venv"
# shellcheck disable=SC1091
source "$HERE/../.venv/bin/activate"
pip install --upgrade pip
pip install -r "$HERE/requirements.txt"

if [ ! -d "$SRC" ]; then
    mkdir -p "$(dirname "$SRC")"
    git clone --recursive --depth 1 --shallow-submodules https://github.com/thliebig/openEMS-Project.git "$SRC"
fi
# The active venv receives the CSXCAD and openEMS Python bindings. FieldForge has its own GUI, so skip AppCSXCAD.
cd "$SRC"
./update_openEMS.sh "$PREFIX" --python --disable-GUI --skip-dep-check

echo "Done. Double-click hfss/FieldForge.command to start FieldForge."
