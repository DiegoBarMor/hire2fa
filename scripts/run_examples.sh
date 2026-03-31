#!/bin/bash
set -euo pipefail

### Run example conversions
if [ ! -d "scripts" ]; then
    echo "Error: script must be run in the repository root directory (where setup.py is located)."
    exit 1
fi

hire2fa testdata/examples/1akx.cg.pdb testdata/examples/1akx.hire2fa.pdb
