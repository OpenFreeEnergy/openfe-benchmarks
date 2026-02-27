# OpenFE Benchmarks Results Submission Format

Purpose
-------
This document defines the structure, required files, metadata schema, and validation expectations for a "submitted result" in this repository. The goal is to make each submission a self‑contained, machine-readable dataset that an external user or automated tool can ingest, validate and cite.

See [OpenFE CLI tutorial](https://docs.openfree.energy/en/latest/tutorials/rbfe_cli_tutorial.html) for more information in how to prepare files for a submission with the plan/gather/run lifecycle and example `gather` output. For visualization and example plots that consume `dg`/`ddg` reports see the [OpenFE plotting tutorial](https://docs.openfree.energy/en/latest/tutorials/plotting.html).

Summary
--------
Each submission is a directory under `results/` that MUST contain: a human- and machine-readable metadata YAML (`submission.yaml`), the tabular outputs produced by `openfe gather` (at minimum `dg` and `ddg` TSVs), and a pointer (DOI or URL) to the full simulation archive (Alchemical Archive / Zenodo) containing the raw JSON/replicate/trajectory assets.

Minimum Required Layout
------------------------
results/<submission-id>/
- `submission.yaml`          # (Required) Metadata + links (see schema)
- `dg.tsv` or `dg.csv`       # (Required) Output of `openfe gather --report dg` (tab- or comma-delimited)
- `ddg.tsv` or `ddg.csv`     # (Required) Output of `openfe gather --report ddg` (tab- or comma-delimited)
- `raw.tsv` or `raw.csv`     # (Recommended) `openfe gather --report raw` (per-leg values)

Note: heavy/binary simulation artifacts (repeat folders, JSON results, trajectories) SHOULD be deposited in a long-term archive (Zenodo/Alchemical Archive) and referenced from `submission.yaml` rather than included in the git repo.

submission.yaml — annotated example (required fields + guidance)
----------------------------------------------------------------
The block below is a single, copy‑pasteable example that documents required and recommended keys inline. Use the exact key names shown; for report filenames the repository accepts `.tsv` or `.csv` (tab- or comma-delimited respectively).

```yaml
# REQUIRED: unique, kebab-case identifier for this submission
submission_id: example-tyk2-2025

# REQUIRED: short descriptive title
title: Tyk2 RBFE benchmark — OpenFE example submission

# REQUIRED: list of contributing authors (name, affiliation; ORCID optional)
authors:
  - name: My Name
    affiliation: Example Lab
    orcid: '0000-0002-XXXX-XXXX'

# REQUIRED: publication/submission date (ISO 8601)
date: 2025-02-10

# REQUIRED: OpenFE version used to produce the gathered reports
openfe_version: 0.8.3

# Recommended but useful: force field and partial charge descriptor
forcefield: openff-2.3.0
partial_charges: am1bcc

# REQUIRED: exact command used to create the reports (provenance)
# can output TSV or CSV; keep the extension consistent with produced files
gather_command: "openfe gather results/ --report dg --report ddg -o final_results.tsv"  # or final_results.csv

# REQUIRED: mapping of logical report names → filenames
# filenames MUST use .tsv or .csv (delimiter must match extension)
reports:
  dg: dg.tsv        # required: absolute/relative path to DG report (dg.csv allowed)
  ddg: ddg.tsv      # required: relative path to DDG report (ddg.csv allowed)
  raw: raw.tsv      # recommended: per-leg values (raw.csv allowed)

# REQUIRED: long-term archive pointer (at least doi or url)
archive:
  doi: 10.5281/zenodo.1234567
  archive_provider: zenodo

# REQUIRED: license for the submission (e.g. CC-BY-4.0)
license: CC-BY-4.0

# RECOMMENDED / OPTIONAL metadata
summary: "RBFE benchmark for TYK2 comparing protocols X and Y; includes dg and ddg summary tables."
keywords: [tyk2, rbfe, benchmark, openfe]
software:
  openmm: 8.0
  openfe: 0.8.3
protocol:
  n_protocol_repeats: 3
  integrator: langevin
  timestep_fs: 2

# notes: any additional provenance or caveats
notes: "Replica counts reduced for testing; full archive available at DOI above."
```

Expected TSV/CSV Contents (What `openfe gather` Writes)
--------------------------------------------------------
Files must be tab- or comma-delimited and use the `.tsv` or `.csv` extension; column headings must match the examples below.

- `dg` report (TSV/CSV): columns typically include at least
  - `ligand` — ligand identifier
  - `DG(MLE) (kcal/mol)` — MLE absolute free energy
  - `uncertainty (kcal/mol)` — reported uncertainty

- `ddg` report (TSV/CSV): columns typically include at least
  - `ligand_i`, `ligand_j`
  - `DDG(i->j) (kcal/mol)` — relative free energy
  - `uncertainty (kcal/mol)`

- `raw` report (TSV/CSV, recommended): per-leg rows with a `leg` column (`vacuum`/`solvent`/`complex`) and the raw DG per leg.

See OpenFE CLI tutorial for concrete examples and column headings: https://docs.openfree.energy/en/latest/tutorials/rbfe_cli_tutorial.html

References
----------
- OpenFE CLI tutorial (gather): https://docs.openfree.energy/en/latest/tutorials/rbfe_cli_tutorial.html