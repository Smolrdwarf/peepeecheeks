#!/bin/bash
# Launch FieldForge on macOS using the repo's virtual environment. Double-click it in Finder.
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export OPENEMS_INSTALL_PATH="${OPENEMS_INSTALL_PATH:-$HOME/opt/openEMS}"
exec "$HERE/../.venv/bin/python" "$HERE/run.py" "$@"
