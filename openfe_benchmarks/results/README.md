# OpenFE Benchmarks Results Submission Format

Purpose
-------
This document defines the structure, required files, metadata schema, and validation expectations for a "submitted result" in this repository. The goal is to make each submission a self‑contained, machine-readable dataset that an external user or automated tool can ingest, validate and cite.

See [OpenFE CLI tutorial](https://docs.openfree.energy/en/latest/tutorials/rbfe_cli_tutorial.html) for more information in how to prepare files for a submission with the plan/gather/run lifecycle and example `gather` output

Summary
--------
Each submission is a directory under `results/` that MUST contain: a human- and machine-readable metadata YAML (`submission.yaml`), the tabular outputs produced by `openfe gather` (at minimum `dg` and `ddg` TSVs), and a pointer (DOI or URL) to the full simulation archive (Alchemical Archive / Zenodo) containing the raw JSON/replicate/trajectory assets.

Minimum Required Layout
------------------------
results/<submission-id>/
- `submission.yaml`          # (Required) Metadata + links (see schema)
- `dg.tsv`                   # (Required) Output of `openfe gather --report dg`
- `ddg.tsv`                  # (Required) Output of `openfe gather --report ddg`
- `raw.tsv`                  # (Recommended) `openfe gather --report raw` (per-leg values)
- `checksums.sha256`         # (Recommended) Checksums for tsv/yaml

Note: heavy/binary simulation artifacts (repeat folders, JSON results, trajectories) SHOULD be deposited in a long-term archive (Zenodo/Alchemical Archive) and referenced from `submission.yaml` rather than included in the git repo.

What must be in `submission.yaml` (schema)
-----------------------------------------
The following minimal schema should be present in every `submission.yaml`. Keys marked **required** must be present.

Required fields
- `submission_id` (str) — **required**: unique id (use kebab-case, e.g. `jensen-tyk2-2025`)
- `title` (str) — **required**: short descriptive title
- `authors` (list) — **required**: list of {name:, affiliation:, orcid: (optional)}
- `contact` (object) — **required**: {name:, email:}
- `date` (ISO 8601) — **required**: date of submission/publication
- `openfe_version` (str) — **required**: OpenFE version used to produce the `gather` outputs
- `gather_command` (str) — **required**: exact `openfe gather` command used (for provenance)
- `reports` (object) — **required**: mapping of report names to filenames; at minimum:
  - `dg: dg.tsv`
  - `ddg: ddg.tsv`
- `archive` (object) — **required**: long-term archive info with at least one of:
  - `doi` (preferred) or `url`
  - `archive_provider` (e.g. `zenodo`, `alchemical-archive`, `institutional`) 
- `license` (str) — **required**: license for the results (e.g. `CC-BY-4.0`)

Recommended/optional fields
- `summary` (str) — short plain‑English description of what is contained and key findings
- `protocol` (object) — summary of protocol settings (e.g. `n_protocol_repeats`, `integrator`, `timestep`, `nonbonded_cutoff`, etc.)
- `software` (object) — versions of important packages (OpenMM, kartograf, nagl, etc.)
- `keywords` (list)
- `associated_publication` (doi/url)
- `compute_environment` (free text — e.g. cluster name, GPU model)
- `checksums` (object) — map of filename → sha256
- `notes` (str)

Example minimal `submission.yaml`
---------------------------------
```yaml
submission_id: example-tyk2-2025
title: Tyk2 RBFE benchmark — OpenFE example submission
authors:
  - name: My Name
    affiliation: Example Lab
    orcid: '0000-0002-XXXX-XXXX'
contact:
  name: My Name
  email: m.name@example.org
date: 2025-02-10
openfe_version: 0.8.3
forcefield: openff-2.3.0
partial_charges:
gather_command: "openfe gather results/ --report dg --report ddg -o final_results.tsv"
reports:
  dg: dg.tsv
  ddg: ddg.tsv
  raw: raw.tsv
archive:
  doi: 10.5281/zenodo.1234567
  archive_provider: zenodo
license: CC-BY-4.0
summary: "RBFE benchmark for TYK2 comparing protocols X and Y; includes dg and ddg summary tables."
keywords: [tyk2, rbfe, benchmark, openfe]
```

Expected TSV Contents (What `openfe gather` Writes)
--------------------------------------------------
- `dg` report (TSV): columns typically include at least
  - `ligand` — ligand identifier
  - `DG(MLE) (kcal/mol)` — MLE absolute free energy
  - `uncertainty (kcal/mol)` — reported uncertainty

- `ddg` report (TSV): columns typically include at least
  - `ligand_i`, `ligand_j`
  - `DDG(i->j) (kcal/mol)` — relative free energy
  - `uncertainty (kcal/mol)`

- `raw` report (TSV, recommended): per-leg rows with a `leg` column (`vacuum`/`solvent`/`complex`) and the raw DG per leg.

See OpenFE CLI tutorial for concrete examples and column headings: https://docs.openfree.energy/en/latest/tutorials/rbfe_cli_tutorial.html

References
----------
- OpenFE CLI tutorial (gather): https://docs.openfree.energy/en/latest/tutorials/rbfe_cli_tutorial.html