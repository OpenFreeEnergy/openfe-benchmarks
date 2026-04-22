Done by Joshua Horton on 2026-02-12 (subset regenerated 2026-04-02)

All ligands extracted from the FreeSolv database version https://github.com/MobleyLab/FreeSolv/releases/tag/v0.52

## Partial Charges

The reference data was generated using the [generate_freesolv_exp_data.py](../../../data_generation/generate_freesolv_exp_data.py) script using the [conda-lock_linux-64.yml](../../../data_generation/conda-lock_linux-64.yml) environment. 
Charges were generated using the [charge_freesolv.py](../../../data_generation/charge_freesolv.py) script using the [conda-lock_linux-64.yml](../../../data_generation/conda-lock_linux-64.yml) environment. 
Some ligands could not be charged with all methods, the following lists the ligands that could not be charged with each method:

Charging failures:
- `am1bccelf10_oe`: mobley_7176248
- `am1bcc_at`: mobley_9741965

## Subsets

Regenerate all subsets: run `run.sh` in `data_generation/`.

### subset_openff_filtered.json
588 neutral solutes retained from the full FreeSolv v0.52 database, all measured in water. Excluded 54 entries: 47 with elements outside the OpenFF chemical space, 5 with disqualifying SMIRKS patterns, and 2 with undefined stereochemistry. Each entry is a unique solute with a single experimental aqueous hydration free energy measurement.

### subset_openff_small.json
188 unique solutes drawn from `subset_openff_filtered`. Selection seeded the 188 solutes present in both the FreeSolv and MNSol filtered pools (by OpenFF SMILES isomorphism).