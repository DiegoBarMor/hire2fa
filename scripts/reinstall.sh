#!/bin/bash
set -eu

### Re-install the package locally

pip uninstall hire2fa -y || true
pip install .
rm -rf build hire2fa.egg-info
