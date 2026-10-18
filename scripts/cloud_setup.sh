#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
python3 -m pip install -e . -c requirements-validation.txt
python3 -c 'import sys, numpy, scipy, sequential_optics; print("Python:", sys.version); print("NumPy:", numpy.__version__, "SciPy:", scipy.__version__); print("Engine:", sequential_optics.__file__)'
