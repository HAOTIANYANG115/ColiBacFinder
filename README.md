# ColiBacFinder

**Reference-based screening of colicin- and microcin-associated loci in bacterial genome assemblies.**

ColiBacFinder searches nucleotide FASTA assemblies with a frozen protein reference database using TBLASTN. It reports bacteriocin/precursor candidates, immunity candidates, component support and genomic neighborhoods, while retaining alternative matches and unresolved types.

**Input:** one genome or plasmid assembly, or a directory of assemblies.

**Start with:** `00_results_report.html` for a readable report and `00_results_summary.tsv` for a spreadsheet-ready summary.

**Interpretation:** sequence-based candidates, not experimentally confirmed bacteriocin production or activity.

## Install with Conda

The published package **0.1.2** is available from the project's [Anaconda.org channel](https://anaconda.org/tiantianyangyang/colibacfinder):

```bash
conda create -n colibacfinder --override-channels \
  -c conda-forge -c bioconda -c tiantianyangyang \
  colibacfinder=0.1.2 -y
conda activate colibacfinder
colibacfinder --version
colibacfinder db-check
```

Expect `0.1.2` and `"checksum_status": "PASS"`. Fresh-environment installation, database verification and English report regeneration have been confirmed for 0.1.2 on Linux x86_64. A real-sample search was previously confirmed with 0.1.1. The package is distributed through a personal channel; it is not a Bioconda release.

**Version status:** 0.1.2 is available through Conda with English summary reports and expanded CLI help. Existing 0.1.1 installations require an upgrade to receive these changes. See [CHANGELOG.md](CHANGELOG.md).

### Dependencies and platforms

- Python >=3.9.
- NCBI BLAST+ >=2.16,<3 (`makeblastdb` and `tblastn` on `PATH`). Conda resolves these dependencies automatically.
- The frozen reference database is bundled; no separate reference download is required.
- No Prokka, ABRicate, Biopython, FAA or GFF input is required.

Linux x86_64 is the tested platform. The Python package is `noarch`, but execution still requires compatible BLAST binaries. Native Windows and macOS operation have not been validated; Windows users can use a Linux environment such as WSL2, subject to local validation.

### Install from source

For development or source-based installation:

```bash
git clone https://github.com/HAOTIANYANG115/ColiBacFinder.git
cd ColiBacFinder
conda env create -f environment.yml
conda activate colibacfinder
python -m pip install --no-deps .
colibacfinder --version
colibacfinder db-check
```

Record the commit with `git rev-parse HEAD` for reproducibility. GitHub access is required if the repository is private.

## Quick start

One assembly:

```bash
colibacfinder run -i genome.fasta -o results_sample -t 8
```

A directory of assemblies (one file per sample):

```bash
colibacfinder run -i genomes/ -o results_batch -t 8
```

The input must contain **nucleotide sequences**, not protein FASTA. Files ending in `.fa`, `.fna` or `.fasta`, optionally `.gz`, are read from the top level of the directory. Other files and subdirectories are ignored. Contigs in one file belong to one sample; do not combine different isolates into a single FASTA.

The output directory must not already exist. Samples are processed sequentially; `-t` controls threads per sample. Inputs are not modified. Interrupted runs are not resumed automatically; use a new output directory for a rerun.

```bash
colibacfinder -h
colibacfinder run -h
```

For long jobs on Linux:

```bash
nohup colibacfinder run -i genomes/ -o results_batch -t 8 > colibacfinder.log 2>&1 &
tail -f colibacfinder.log
```

Exit codes: `0`, all samples completed; `1`, at least one sample failed; `2`, argument or preparation error; `130`, interrupted.

## Read the results

Open **`00_results_report.html`** in a browser, or import **`00_results_summary.tsv`** as a tab-delimited table. Each row represents one sample.

| Summary field | Meaning |
|---|---|
| `run_status` | Whether the analysis completed; failed samples are not negatives |
| `result` | Near-reference bacteriocin/precursor, immunity only, weak/partial candidates only, or no primary candidate |
| `main_candidates` | Systems supported by a near-reference bacteriocin or precursor hit; unresolved alternatives are retained |
| `component_support` | Number of supported components out of those included in the current model |
| `immunity_candidates` | Near-reference immunity labels, which may include multiple labels for one locus |
| `weak_or_partial_loci` | Weak/partial system anchor loci; inspect detailed tables for all candidate regions |
| `notes` | Ambiguity, warnings or errors |

For example, `Colicin_M: 2/2 [bacteriocin_immunity_only]` means that both modeled components (bacteriocin and immunity) have support. It **does not** establish an intact, expressed or active system. Names separated by `/` are competing types at the same locus, not additional genes. An immunity hit alone does not establish bacteriocin production.

A negative screening result means **no qualifying candidate in this database at these thresholds**, not proven biological absence.

### Detailed evidence

| Output | Contents |
|---|---|
| `01_sample_summary.tsv` | Technical per-sample counts, labels and failures |
| `02_all_HSPs.tsv` | Retained search hits, classifications, coordinates and reference annotations |
| `03_candidate_loci.tsv` | Merged candidate loci, best/alternative matches, score differences and edge flags |
| `04_locus_evidence.tsv` | Individual matches supporting each locus |
| `05_proximity_pairs.tsv` | Bacteriocin/immunity proximity, overlap and strand; proximity is not functional confirmation |
| `06_system_candidates.tsv` | System candidates and component-support status |
| `07_system_components.tsv` | Evidence for individual components |
| `08_label_frequency.tsv` | Sample-level label frequencies among completed samples; categories can overlap |
| `samples/S00001/` and subsequent directories | Raw hits, logs, contig ID maps and candidate neighborhood sequences |
| `input_manifest.tsv` | Sample IDs, input paths, checksums and work directories |
| `parameters.json`, `system_models.json` | Versions, screening rules and component-model definitions |
| `reference_snapshot/` | Reference sequences, source accessions and evidence metadata used in the run |

Labels and HSPs are not independent gene counts. Candidate loci are not confirmed full-length proteins; coding integrity is not assessed.

## Screening rules

The search retains hits at E-value <=1e-3. Screening requires E-value <=1e-5 and applies the following classes in order:

| Class | Identity | Reference span | Minimum alignment length |
|---|---:|---:|---:|
| `near_reference` | >=80% | >=80% | No additional minimum |
| `remote` | >=35% | >=70% | 40 aa |
| `partial` | >=35% | >=10% | 60 aa |

Reference span is the percentage of reference coordinates covered by **one HSP**, not cumulative coverage across HSPs. Thresholds are exploratory and have not been calibrated on an independent biological truth set. A high-scoring shared domain alone does not establish a complete bacteriocin or a precise subtype.

Useful options:

| Option | Default | Purpose |
|---|---:|---|
| `--identity`, `--coverage` | 80, 80 | Near-reference screening thresholds (%) |
| `--evalue` | 1e-5 | Screening E-value |
| `--search-evalue` | 1e-3 | Search E-value |
| `--remote-identity`, `--remote-coverage` | 35, 70 | Remote-candidate thresholds (%) |
| `--partial-identity`, `--partial-coverage` | 35, 10 | Partial-candidate thresholds (%) |
| `--neighbor-bp` | 10000 | Maximum component-to-anchor gap |
| `--flank-bp` | 5000 | Sequence extracted on each side of a candidate |
| `--db` | Bundled release | External compatible database directory |
| `--circular-map` | None | TSV with `sample` and `contig` columns for explicitly closed molecules |
| `--keep-db` | Off | Retain per-sample BLAST databases |

Overlapping precursor/immunity hits, including those in Microcin V, are retained with explicit overlap lengths. Shared export or maturation proteins do not count as bacteriocin precursors. Model scopes differ among systems; inspect `system_models.json` before interpreting component completeness.

Circular declarations enable origin-aware neighborhood distances and extraction. They do not join TBLASTN HSPs across the origin. Declare circularity only when closure is supported by assembly evidence.

## Regenerate reports

Generate the English summary from an existing unified-version run after installing version 0.1.2:

```bash
colibacfinder summarize --run results_batch
```

This replaces only the two `00_` reports, preserves the detailed evidence and uses the existing classifications. No BLAST rerun is needed. Version 0.1.1 also supports this command but produces its original summary language.

To reclassify retained raw hits into a **new** report directory:

```bash
colibacfinder report --from-run results_batch -o results_reclassified
```

This checks available input checksums and database fingerprints. Threshold options can be supplied; see `colibacfinder report -h`. Legacy batch-script runs lack some provenance information and are flagged accordingly. Counts need not match those of the earlier exploratory scripts.

## Reference database and reproducibility

Bundled release: `bacteriocin_unified_v0.1_20261020_194442_905507`.

- 84 distinct protein sequences and 92 annotation relationships.
- Main, candidate, fragment and auxiliary references retain their evidence tiers.
- 39 distinct source CDS sequences, linked to 36 proteins; protein sequences without CDS links are not reverse-translated.
- Original reference bytes and SHA256 checksums are preserved across this presentation update.

Accessions and evidence links are available in `reference_snapshot/metadata/reference_links.tsv`; detailed provenance is in `reference_snapshot/source_metadata/`. Historical source metadata are preserved as supplied, including original-language provenance text. Report the software version, database fingerprint and thresholds in research outputs.

## Validation and limitations

The published 0.1.1 package passed fresh-environment installation, database verification and a real plasmid run reported by the maintainer. The published 0.1.2 package subsequently passed fresh-environment installation, checksum verification and English summary regeneration from that retained run. No new BLAST search was performed in this 0.1.2 check. These are deployment and reporting checks, not independent accuracy validation. See [docs/VALIDATION.md](docs/VALIDATION.md) for software checks and remaining evaluation needs.

The software does not identify species, determine plasmid versus chromosome location, verify intact ORFs, establish expression/activity or guarantee exact subtype calls. Fragmented assemblies, short peptides, divergent sequences and reference gaps can cause missed or ambiguous hits. The search is protein-to-translated-DNA; there is no separate nucleotide BLAST detection path.

For independent assessment:

```bash
colibacfinder evaluate --run results_batch --truth truth.tsv -o evaluation
```

The tab-delimited truth table requires `sample`, `target`, `truth` (0/1) and `evidence`. Supported targets are `ANY_CORE`, `ANY_NONEXPORT_IMMUNITY`, `system:Microcin_V`, `system:Microcin_J25` and `system:Microcin_B17`. Use independent evidence for the same endpoint being evaluated. Do not label unknowns as negative or interpret a negative inhibition assay as proof of gene absence. Failed/missing samples are excluded explicitly. Reference self-matches and synthetic controls are software checks, not estimates of biological sensitivity or specificity.

## Development and support

- [Build and publish a Conda package](docs/CONDA_RELEASE.md)
- [Validation record](docs/VALIDATION.md)
- [Version history](CHANGELOG.md)
- [Report an issue](https://github.com/HAOTIANYANG115/ColiBacFinder/issues): include the software version, database-check output, command and relevant logs; do not include credentials or confidential genomic data.

No software license has yet been specified in this repository. Third-party reference records retain their original terms; package availability alone does not grant redistribution rights. No peer-reviewed software citation is currently specified; identify the repository, version and database release when describing use.
