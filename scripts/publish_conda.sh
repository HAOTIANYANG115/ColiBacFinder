#!/usr/bin/env bash
set -euo pipefail
if [[ $# != 2 ]]; then
  echo "Usage: bash scripts/publish_conda.sh ANACONDA_USERNAME PACKAGE_PATH" >&2
  exit 2
fi
CHANNEL=$1
PACKAGE=$2
[[ -f "$PACKAGE" ]] || { echo "Package not found" >&2; exit 2; }
# Authenticate locally with anaconda login; never place tokens in this script.
conda build --test "$PACKAGE" --override-channels -c conda-forge -c bioconda
anaconda upload --user "$CHANNEL" "$PACKAGE"
printf 'After upload, users can run: conda create -n colibacfinder --override-channels -c %s -c conda-forge -c bioconda colibacfinder=0.1.1\n' "$CHANNEL"
