# Conda distribution and release guide

## Published version

Version **0.1.2** is available at https://anaconda.org/tiantianyangyang/colibacfinder.
Fresh installation, database integrity and English report regeneration were verified on Linux x86_64. This is a personal channel, not a Bioconda release.

```bash
conda create -n colibacfinder --override-channels \
  -c conda-forge -c bioconda -c tiantianyangyang \
  colibacfinder=0.1.2 -y
conda activate colibacfinder
colibacfinder --version
colibacfinder db-check
colibacfinder run -i genomes/ -o new_results -t 8
```

Python >=3.9 and BLAST >=2.16,<3 are runtime dependencies. Users do not need Prokka, ABRicate, Biopython or the build tools. The Python package is noarch, but supported execution depends on BLAST availability; Linux x86_64 is the tested platform.

## Build a package (maintainers, Linux)

```bash
conda create -n cbf-build --override-channels -c conda-forge \
  python=3.11 conda-build anaconda-client=1.14.1 pip setuptools wheel -y
conda activate cbf-build
cd /path/to/ColiBacFinder
bash scripts/check_package.sh
bash scripts/build_conda.sh
```

The scripts do not upload automatically. Conda tests CLI entry points, database integrity and BLAST dependencies in its package test environment. Artifacts are written to `conda-dist/noarch/`. Use the actual generated filename; the current version is 0.1.2; future releases require synchronized version updates. Stop if tests fail.

If the system `/tmp` is read-only, select a writable directory before building:

```bash
mkdir -p "$HOME/.cache/colibacfinder/tmp"
export TMPDIR="$HOME/.cache/colibacfinder/tmp"
export TMP="$TMPDIR"
export TEMP="$TMPDIR"
```

Never regenerate a checksum manifest merely to hide a packaging mismatch. The frozen reference files must preserve their original bytes; `.gitattributes` disables text conversion for the database directory.

## Authenticate and upload

For legacy anaconda-client 1.14.1, the following environment variable selects the standalone client when newer CLI plugins are also installed:

```bash
export ANACONDA_CLIENT_FORCE_STANDALONE=1
anaconda login
anaconda whoami
bash scripts/publish_conda.sh tiantianyangyang /absolute/path/to/colibacfinder-0.1.2-py_0.conda
```

Use Anaconda.org credentials, not the password of a social-login provider. If token authentication is used, obtain an Anaconda.org token with the necessary repository/API permissions and supply it locally through the client's supported authentication mechanism. Never commit or share passwords or tokens.

The upload script tests the artifact first and does not force replacement of an existing package. Confirm the intended account, package version and distribution permissions before publication. This repository does not yet specify a software license; preserve third-party reference terms.

## Verify the published artifact

To reproduce the installation check, install from the channel in a new environment:

```bash
conda create -n cbf-release-test-012 --override-channels \
  -c conda-forge -c bioconda -c tiantianyangyang \
  colibacfinder=0.1.2 -y
conda activate cbf-release-test-012
colibacfinder --version
colibacfinder db-check
colibacfinder run -i /path/to/test_genome.fasta -o /path/to/new_test_results -t 8
```

Check the expected version, checksum PASS, successful run status and English summary output. Update the README's published version only after upload and installation verification. A GitHub source update alone does not change an existing Conda package.
