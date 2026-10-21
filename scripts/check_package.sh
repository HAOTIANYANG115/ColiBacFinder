#!/usr/bin/env bash
set -euo pipefail
PROJECT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
CHECK=$(mktemp -d "${TMPDIR:-/tmp}/cbf-package.XXXXXXXX")
trap 'rm -rf "$CHECK"' EXIT
# A separate source copy prevents stale build/ artifacts from entering the wheel.
python - "$PROJECT" "$CHECK/source" <<'PY'
import shutil,sys
shutil.copytree(sys.argv[1],sys.argv[2],ignore=shutil.ignore_patterns('.git','build','dist','*.egg-info','__pycache__','conda-dist'))
PY
python -m pip wheel --no-deps --no-build-isolation "$CHECK/source" -w "$CHECK/wheels"
python -m venv "$CHECK/env"
"$CHECK/env/bin/python" -m pip install --no-deps "$CHECK"/wheels/colibacfinder-*.whl
cd "$CHECK"
unset PYTHONPATH
"$CHECK/env/bin/colibacfinder" --version
"$CHECK/env/bin/colibacfinder" db-check
"$CHECK/env/bin/colibacfinder" summarize -h
printf 'PASS: isolated wheel installation and original database checksums\n'
