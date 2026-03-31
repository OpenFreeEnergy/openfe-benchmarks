Done by Joshua Horton on 2026-02-12

All ligands extracted from the FreeSolv database version https://github.com/MobleyLab/FreeSolv/releases/tag/v0.52

## Subsets

Of the 603 systems imported from the FreeSolv database, two subsets were made:

- Filtered systems: 588
- Overlapping FreeSolv solutes with MNSol: 189
- Final aligned subset rows: 201
- Missing chemical environments: 18

Missing environments:

CarbonylHydrate, Hemiaminal, Aminal, Thioacetal, Cyanate, Isocyanate, AlkylFluoride, ArylFluoride, Thioaldehyde, Thioketone, Thiourea, Thiocyanate, Isothiocyanate, ThiocarboxylicAcid, ThiocarboxylicAcidEster, SulfonicAcidEster, PhosphonicAcid, PhosphoricAcid

Subset files and regeneration script:

- subset_openff_mnsolv_aligned_filtered.json
- subset_openff_mnsolv_aligned_small.json
- ../../../data_generation/define_freesolv_mnsolv_openff_subsets.py

## Notes

The reference data was generated using the [generate_freesolv_exp_data.py](../../../data_generation/generate_freesolv_exp_data.py) script using the [conda-lock_linux-64.yml](../../../data_generation/conda-lock_linux-64.yml) environment. 
Charges were generated using the [charge_freesolv.py](../../../data_generation/charge_freesolv.py) script using the [conda-lock_linux-64.yml](../../../data_generation/conda-lock_linux-64.yml) environment. 
Some ligands could not be charged with all methods, the following lists the ligands that could not be charged with each method:

- am1bccelf10_oe:
    - mobley_7176248

- am1bcc_at:
    - mobley_9741965