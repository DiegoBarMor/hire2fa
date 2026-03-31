#!/bin/bash
set -euo pipefail

### Build and upload the package to PyPI
if [ ! -d "scripts" ]; then
    echo "Error: script must be run in the repository root directory (where setup.py is located)."
    exit 1
fi

# pip install build twine
python3 -m build
twine upload dist/*
rm -rf build dist hire2fa.egg-info
