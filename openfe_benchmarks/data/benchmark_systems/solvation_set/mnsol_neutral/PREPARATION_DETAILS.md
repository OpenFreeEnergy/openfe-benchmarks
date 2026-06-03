Done by Jennifer Clark on 2026-02-25 (subset regenerated 2026-05-21 by Jennifer Clark)

Systems extracted from the MNSol database (https://doi.org/10.13020/3eks-j059).

- Licensing prevented including MNSol experimental values. After properly obtaining the dataset, a user can locally generate the needed experimental data file with the [generate_mnsol_data.py](../../../data_generation/generate_mnsol_data.py) script.
- Experimental uncertainties are set to 0.2 kcal/mol for all neutral entries, following the recommendation in the MNSol documentation.

## Filtering

The reference data was generated using the [generate_mnsol_data.py](../../../data_generation/generate_mnsol_data.py) script. Entries were excluded if:

- The solute or solvent name was not present in `mnsol-name-to-smiles.json`
- Molecules with ambiguous isomeric structure were excluded including: 'bromotoluene', 'chlorotoluene', 'dichloroethane', 'fluoroctane', 'trimethylbenzene'
- The solute is a radical (name contains "radical")
- The solute is a dimer (MNSol Level1=14)
- The net charge is non-zero
- The solvent has undefined stereochemistry

Conformer generation failed for: 6,7,8,9,10,10-hexachloro-1,5,5a,6,9,9a-hexahydro-6,9-methano-2,4,3-benzodioxathiepine-3-oxide(endosulfanalpha)

After filtering, **2467 systems** are retained.

## Partial Charges

Charges were generated using the [charge_mnsol.py](../../../data_generation/charge_mnsol.py) script using the [conda-lock_linux-64.yml](../../../data_generation/conda-lock_linux-64.yml) environment. 
Some ligands could not be charged with all methods, the following lists the ligands that could not be charged with each method:

- `am1bcc_oe`, `am1bcc_at`, `am1bccelf10_oe`: hydrogen
- `nagl_off`: hydrogen, tetramethylsilane

## Subsets

Regenerate all subsets: run `python define_freesolv_mnsol_openff_subsets.py` in `data_generation/`.

### subset_openff_filtered.json
These subsets are independent of subsets used for Sage < 2.4.

1570 systems retained from ~3037 MNSol database entries, covering 236 unique solutes across multiple solvent environments. Excluded: 402 entries involving water; 71 self-solvation entries (solute = solvent); 47 with solutes containing out-of-scope elements; 40 with solvents containing out-of-scope elements; 3 entries with solutes bearing disqualifying SMIRKS; 286 with solvents bearing disqualifying SMIRKS; and 1 with undefined solvent stereochemistry. Entries may satisfy multiple exclusion criteria; totals are not additive.

### subset_openff_small.json
These subsets are independent of subsets used for Sage < 2.4.

400 systems, 195 unique solutes, with a per-solute cap of 3. Solute repetition: 79 appear once (62 pool-exhausted, 17 selection-limited), 27 appear twice (19 pool-exhausted, 8 selection-limited), and 89 appear three times (17 pool-exhausted, 72 at cap). Selection seeded the 188 MNSol–FreeSolv overlap solutes (Phase 1, filling the full 400-system budget), then applied Tanimoto pair-space rebalancing with 500 swaps (486 same-solute, 14 cross-solute). ChemicalEnvironment coverage fill (Phase 2) was not reached as the budget was consumed by Phase 1. 187 / 195 MNSol solutes (95.9%) overlap with the FreeSolv filtered pool. 24 ChemicalEnvironments have fewer than 3 representatives: 19 absent from the filtered pool (0 representatives), 1 pool-limited (Disulfide: all 2 available representatives selected), and 4 budget-limited (CarboxylicAcidSecondaryAmide, Sulfone, Sulfoxide, TertiaryAmine: ≥3 in pool, not all selected).