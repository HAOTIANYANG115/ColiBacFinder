# Changelog

## 0.1.2 — published to Anaconda.org on 2026-10-22

- Use English for summary TSV conclusions, HTML reports and interpretation notes.
- Expand command-line help with input requirements, parameter descriptions and defaults.
- Replace the README and release guide with English documentation for international researchers.
- Document the published 0.1.1 Conda channel and its successful Linux installation and real-sample check.
- Preserve all screening thresholds, detailed table schemas, reference bytes and evidence tiers.
- Verify fresh Conda installation, database checksum PASS and English summary regeneration on Linux.
- Keep 0.1.1 available; the 0.1.2 changes do not update already installed packages automatically.

## 0.1.1 — published to Anaconda.org on 2026-10-22

- Automatically generate `00_results_summary.tsv` and `00_results_report.html`.
- Add `summarize --run` to generate summaries without rerunning BLAST.
- Keep weak/partial candidates separate from main candidates and retain competing types at one anchor.
- Distinguish failed analyses from no-detection results.
- Preserve original database bytes and checksums; add isolated wheel installation checks.
- Add a Conda recipe, build/upload scripts and explicit BLAST runtime dependencies.
- Retain the original screening thresholds and reference evidence tiers.
- Publish `colibacfinder-0.1.1-py_0.conda` to the `tiantianyangyang` channel.
