Done by Joshua Horton on 2026-02-12 (subset regenerated 2026-05-21 Jennifer Clark)

All ligands extracted from the FreeSolv database version https://github.com/MobleyLab/FreeSolv/releases/tag/v0.52

## Partial Charges

The reference data was generated using the [generate_freesolv_exp_data.py](../../../data_generation/generate_freesolv_exp_data.py) script using the [conda-lock_linux-64.yml](../../../data_generation/conda-lock_linux-64.yml) environment. 
Charges were generated using the [charge_freesolv.py](../../../data_generation/charge_freesolv.py) script using the [conda-lock_linux-64.yml](../../../data_generation/conda-lock_linux-64.yml) environment. 
Some ligands could not be charged with all methods, the following lists the ligands that could not be charged with each method:

Charging failures:
- `am1bccelf10_oe`: mobley_7176248
- `am1bcc_at`: mobley_9741965

## Subsets

Regenerate all subsets: run `python define_freesolv_mnsol_openff_subsets.py` in `data_generation/`.

### subset_openff_filtered.json
These subsets are independent of subsets used for Sage < 2.4.

588 neutral solutes retained from the full FreeSolv v0.52 database, all measured in water. Excluded 54 entries: 47 with elements outside the OpenFF chemical space, 5 with disqualifying SMIRKS patterns, and 2 with undefined stereochemistry. Each entry is a unique solute with a single experimental aqueous hydration free energy measurement.

### subset_openff_small.json
206 unique solutes drawn from `subset_openff_filtered`. Selection seeded the 188 solutes present in both the FreeSolv and MNSol filtered pools (by OpenFF SMILES isomorphism), then added 18 further solutes via ChemicalEnvironment coverage fill (Tanimoto MaxMin, solute-only environments). 188 / 206 FreeSolv solutes (91.3%) overlap with the MNSol filtered pool. 23 ChemicalEnvironments are not represented: 15 absent from the filtered pool, 8 pool-limited (fewer than 3 representatives available).