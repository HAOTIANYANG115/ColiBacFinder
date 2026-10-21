#!/usr/bin/env bash
set -euo pipefail
PROJECT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$PROJECT"
command -v conda >/dev/null
python -m colibacfinder db-check
python -m unittest discover -s tests -p 'test_*.py'
# Always test the built package; never upload automatically.
conda build recipe --override-channels -c conda-forge -c bioconda --no-anaconda-upload --output-folder "$PROJECT/conda-dist"
printf '\nConda build and package tests passed. Packages in: %s/conda-dist/noarch\n' "$PROJECT"
