# Validation record — 2026-10-21

## Completed

- 25 standard-library unit tests passed. Includes reference integrity, input validation, threshold boundaries, locus ambiguity/deduplication, V overlap, shared component safeguards, fragment handling, circular distance/extraction and failed-case logic. See `unit_tests.log` and `tests/test_rules.py`.
- NCBI BLAST+ 2.16.0 actual integration test passed. Built a nucleotide database and searched with the packaged 84-protein reference. A synthetic construct assembled from reference CDS gave the expected V component-model call; a simple repeat sequence gave no near-reference core call. Invalid FASTA was marked FAILED and excluded from frequency and evaluation denominators.
- Re-reporting the retained hits reproduced candidate locus output byte-for-byte. Existing output directories were protected against overwrite. See `blast_smoke.log` and `tests/smoke_blast.py`.
- Wheel installed in an isolated virtual environment and bundled reference checksum verification passed.

These tests establish exercised software behavior only. Reference-derived constructs are not independent positive controls, and repeat sequences are not representative biological negative controls. Their evaluation metrics must not be used as sensitivity/specificity estimates.

## Not completed

- Independent biological truth-set evaluation, threshold calibration or taxon-wide specificity assessment.
- Re-running the user's 524, 96, 11 or APEC genome cohorts with this release; the supplied source archive contains scripts/database, not those genomic inputs.
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
