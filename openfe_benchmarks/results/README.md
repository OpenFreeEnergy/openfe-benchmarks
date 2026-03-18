# OpenFE Benchmarks Results Submission Format

Purpose
-------
This document defines the structure, required files, metadata schema, and validation expectations for a "submitted result" in this repository. The goal is to make each submission a self‑contained, machine-readable dataset that an external user or automated tool can ingest, validate and cite.

~~See [OpenFE CLI tutorial](https://docs.openfree.energy/en/latest/tutorials/rbfe_cli_tutorial.html) for more information in how to prepare files for a submission with the plan/gather/run lifecycle and example `gather` output. For visualization and example plots that consume `dg`/`ddg` reports see the [OpenFE plotting tutorial](https://docs.openfree.energy/en/latest/tutorials/plotting.html).~~

Summary
--------
Each submission is a directory under `results/` that MUST contain: a human- and machine-readable metadata YAML (`submission.yaml`), a `computational_results.json` containing `dg` and `ddg` information with units for each benchmark system, and a pointer (DOI or URL) to the full simulation archive (Alchemical Archive / Zenodo) containing the raw JSON/replicate/trajectory assets.

Minimum Required Layout
------------------------
results/<submission-id>/
- `submission.yaml`          # (Required) Metadata + links (see schema)
- `computational_results.json`    # (Required) Ligand, dg, ddg, and minimal mbar information

Note: heavy/binary simulation artifacts (repeat folders, openfe JSON results, trajectories, AlchemicalArchives) SHOULD be deposited in a long-term archive (Zenodo/Alchemical Archive) and referenced from `submission.yaml` rather than included in the git repo.

submission.yaml — annotated example (required fields + guidance)
----------------------------------------------------------------
The block below is a single, copy‑pasteable example that documents required and recommended keys inline. Use the exact key names shown.

```yaml
# REQUIRED: unique, kebab-case identifier for this submission
submission_id: example-tyk2-2025

# REQUIRED: short descriptive title
title: Tyk2 RBFE benchmark — OpenFE example submission

# REQUIRED: comprehensive description of what networks the submission calculated and what is / is not comparable to the norm.
description: This submission presents a RBFE benchmark study for JACS TYK2. The calculations employ the standard OpenFE protocol (SOMD with AMBER-based forcefield, AM1-BCC charges) using 3 independent repeats. Note: this submission focuses on relative free energies only; absolute binding free energies (DG) are derived from a reference ligand and should not be interpreted as stand-alone predictions.

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

# REQUIRED: mapping of report names
reports: computational_results.json

# REQUIRED: long-term archive pointer (at least doi or url)
archive:
  doi: 10.5281/zenodo.1234567
  archive_provider: zenodo

# REQUIRED: license for the submission (e.g. CC-BY-4.0)
license: CC-BY-4.0

# RECOMMENDED / OPTIONAL metadata
summary: "RBFE benchmark for TYK2 comparing protocols X and Y; includes dg and ddg summary tables."
keywords: [tyk2, rbfe, benchmark, openfe]
```

Expected `computational_results.json`
-----------------------------------------

The file is a JSON object with two top-level arrays: `"DG"` and `"DDG"`.

### Unit objects

Both `DG`/`DG_uncertainty` and `DDG`/`DDG_uncertainty` values are encoded as a unit-aware object:

```json
{
  "magnitude": -0.312,
  "unit": "kilocalories_per_mole",
  ":is_custom:": true,
  "pint_unit_registry": "openff_units"
}
```

| key | purpose |
|-----|---------|
| `magnitude` | numerical value |
| `unit` | must be `"kilocalories_per_mole"` |
| `":is_custom:"` | always `true`; marks this as an openff-units Quantity |
| `pint_unit_registry` | always `"openff_units"`; identifies the unit registry used for deserialisation |

---

### `DG` array

One entry per ligand per system. These are the **per-ligand absolute binding free energies** derived from a network-wide Maximum Likelihood Estimate (MLE) across all DDG edges.

```json
{
  "ligand": "ejm_50",
  "DG": { "magnitude": 0.291, "unit": "kilocalories_per_mole", ":is_custom:": true, "pint_unit_registry": "openff_units" },
  "DG_uncertainty": { "magnitude": 0.099, "unit": "kilocalories_per_mole", ":is_custom:": true, "pint_unit_registry": "openff_units" },
  "system_group": "jacs_set",
  "system_name": "tyk2",
  "source": "MLE"
}
```

| key | purpose |
|-------|---------|
| `ligand` | ligand identifier matching the SDF/benchmark data |
| `DG` | estimated absolute binding free energy (kcal/mol); note this is an MLE-shifted relative quantity, not a true absolute binding free energy — its absolute value is only meaningful relative to the reference ligand within the same network |
| `DG_uncertainty` | 1-sigma propagated uncertainty from the MLE fit |
| `system_group` | benchmark set the system belongs to (e.g. `jacs_set`, `charge_annihilation_set`) |
| `system_name` | protein target (e.g. `tyk2`, `p38`) |
| `source` | estimation method; currently always `"MLE"` |

---

### `DDG` array

One entry per perturbation edge. These are the **raw pairwise relative free energies** computed directly from simulations, before any network MLE.

```json
{
  "ligand_a": "ejm_50",
  "ligand_b": "ejm_42",
  "system_group": "jacs_set",
  "system_name": "tyk2",
  "repeats": 3,
  "DDG": { "magnitude": -0.312, "unit": "kilocalories_per_mole", ":is_custom:": true, "pint_unit_registry": "openff_units" },
  "DDG_uncertainty": { "magnitude": 0.185, "unit": "kilocalories_per_mole", ":is_custom:": true, "pint_unit_registry": "openff_units" },
  "DGs_complex": [-19.161, -19.111, -18.822],
  "DGs_solvent": [-18.572, -18.830, -18.756],
  "Complex_smallest_mbar_overlaps": [0.134, 0.131, 0.130],
  "Complex_smallest_replica_mixing": [0.089, 0.096, 0.092],
  "Solvent_smallest_mbar_overlaps": [0.138, 0.138, 0.138],
  "Solvent_smallest_replica_mixing": [0.109, 0.105, 0.110]
}
```

| key | purpose |
|-------|---------|
| `ligand_a` / `ligand_b` | the two ligands defining the perturbation edge; DDG = DG(ligand_b) − DG(ligand_a) |
| `system_group` / `system_name` | same meaning as in `DG` entries |
| `repeats` | number of independent simulation repeats; determines the length of each per-repeat array below |
| `DDG` | mean relative binding free energy across repeats (kcal/mol) |
| `DDG_uncertainty` | standard error of the mean across repeats |
| `DGs_complex` | per-repeat ΔG of the alchemical transformation in the **protein–ligand complex** (kcal/mol); length equals `repeats` |
| `DGs_solvent` | per-repeat ΔG of the alchemical transformation in the **solvent** leg (kcal/mol); length equals `repeats` |
| `Complex_smallest_mbar_overlaps` | per-repeat **minimum pairwise MBAR overlap** across all alchemical windows in the complex leg; values near 0 indicate a window spacing that is too coarse and may compromise free energy accuracy |
| `Complex_smallest_replica_mixing` | per-repeat **minimum replica-exchange acceptance rate** between neighbouring windows in the complex leg; values near 0 indicate poor mixing and potentially unconverged sampling |
| `Solvent_smallest_mbar_overlaps` | same as `Complex_smallest_mbar_overlaps` but for the solvent leg |
| `Solvent_smallest_replica_mixing` | same as `Complex_smallest_replica_mixing` but for the solvent leg |


