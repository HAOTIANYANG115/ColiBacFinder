# Validation record

## Completed

- 25 standard-library unit tests passed. Includes reference integrity, input validation, threshold boundaries, locus ambiguity/deduplication, V overlap, shared component safeguards, fragment handling, circular distance/extraction and failed-case logic. See `unit_tests.log` and `tests/test_rules.py`.
- NCBI BLAST+ 2.16.0 actual integration test passed. Built a nucleotide database and searched with the packaged 84-protein reference. A synthetic construct assembled from reference CDS gave the expected V component-model call; a simple repeat sequence gave no near-reference core call. Invalid FASTA was marked FAILED and excluded from frequency and evaluation denominators.
- Re-reporting the retained hits reproduced candidate locus output byte-for-byte. Existing output directories were protected against overwrite. See `blast_smoke.log` and `tests/smoke_blast.py`.
- Wheel installed in an isolated virtual environment and bundled reference checksum verification passed.

These tests establish exercised software behavior only. Reference-derived constructs are not independent positive controls, and repeat sequences are not representative biological negative controls. Their evaluation metrics must not be used as sensitivity/specificity estimates.

## Not completed

- Independent biological truth-set evaluation, threshold calibration or taxon-wide specificity assessment.
- Independent cohort-wide reassessment of the 524, 96, 11 and APEC datasets. Earlier exploratory runs and maintainer-reported checks do not constitute a held-out benchmark.
- Performance/runtime benchmarks on hundreds of complete genomes.
- Experimental activity confirmation, intact ORF assessment, cross-origin HSP joining, or exact colicin subtype disambiguation.

## Independent evaluation protocol

Freeze this package, reference fingerprint and thresholds before testing. Assemble held-out samples with traceable evidence for the exact endpoint: genomic core presence, immunity presence or enumerated system components. Unknowns must remain unknown. Avoid counting near-clonal copies as independent evidence; split by lineage/reference relatedness and retain a separate development set.

Use the `evaluate` command with evidence-linked truth labels. Failed analyses are excluded explicitly, not counted as negatives. Review false positives and false negatives individually and stratify by system, assembly quality, reference tier and sequence relatedness. Functional inhibition and gene presence are different endpoints and require separate validation.

## Reproduce software checks

From the extracted source directory with BLAST+ on PATH:

```bash
python -m unittest discover -s tests -v
PYTHONPATH=. python tests/smoke_blast.py
```

The integration test creates and removes its own temporary directory only; it does not modify user genomes.

## Published Conda 0.1.1: maintainer-reported check, 2026-10-22

The maintainer supplied terminal output showing installation in a fresh `cbf-release-test` Linux environment, version 0.1.1, a working CLI, and database checksum PASS (84 proteins; 92 annotation links). The installed package used Python 3.15. A real-input run on `pe23830-PAI.fa` completed with one near-reference core locus, `Colicin_M`, bacteriocin/immunity support of 2/2 and one weak/partial locus. This reproduced the earlier candidate interpretation. The raw genomic input was not independently inspected during this documentation update.

The published artifact is `colibacfinder-0.1.1-py_0.conda` in the `tiantianyangyang` Anaconda.org channel. These observations support installation and execution reproducibility, not biological accuracy or activity.

## English source update 0.1.2

All 31 unit tests passed, including the 29 existing tests and new checks for formatted CLI help and English HTML with escaped sample names. Database content, thresholds, table field names and system rules are unchanged. HTML and summary TSV presentation text now use English; the HTML language is `en`. CLI help describes parameters and displays defaults. The 0.1.2 Conda artifact was subsequently built and tested successfully; 0.1.1 installations do not acquire these changes automatically.

## Published Conda 0.1.2: maintainer-reported check, 2026-10-22

Terminal output confirmed a successful Conda build and package tests, followed by installation in a fresh `cbf-release-test-012` Linux environment (Python 3.15). The CLI reported 0.1.2 and database checksum PASS, with 84 proteins and 92 annotation links. The manifest and protein FASTA fingerprints were unchanged.

Running `summarize` on the retained 0.1.1 `pe23830-PAI` run generated English reports. The summary retained `Colicin_M`, 2/2 bacteriocin/immunity component support and one weak/partial locus. This verifies installed report generation and preservation of the earlier classifications; it is not a new 0.1.2 BLAST run or independent biological validation.
