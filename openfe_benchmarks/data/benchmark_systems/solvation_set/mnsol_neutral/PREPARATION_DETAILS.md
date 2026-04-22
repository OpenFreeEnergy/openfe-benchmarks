Done by Jennifer Clark on 2026-02-25 (subset regenerated 2026-04-02)

Systems extracted from the MNSol database (https://doi.org/10.13020/3eks-j059).

## Filtering

The reference data was generated using the [generate_mnsol_data.py](../../../data_generation/generate_mnsol_data.py) script. Entries were excluded if:

- The solute or solvent name was not present in `mnsol-name-to-smiles.json`
- Molecules with ambiguous isomeric structure were excluded including: 'bromotoluene', 'chlorotoluene', 'dichloroethane', 'fluoroctane', 'trimethylbenzene'
- The solute is a radical (name contains "radical")
- The solute is a dimer (MNSol Level1=14)
- The net charge is non-zero

Conformer generation failed for: 6,7,8,9,10,10-hexachloro-1,5,5a,6,9,9a-hexahydro-6,9-methano-2,4,3-benzodioxathiepine-3-oxide(endosulfanalpha)

After filtering, **2467 systems** are retained.

## Partial Charges

- Licensing prevented including MNSol experimental values. After properly obtaining the dataset, a user can locally generate the needed experimental data file with the [generate_mnsol_data.py](../../../data_generation/generate_mnsol_data.py) script.
- Experimental uncertainties are set to 0.2 kcal/mol for all neutral entries, following the recommendation in the MNSol documentation.

Charges were generated using the [charge_mnsol.py](../../../data_generation/charge_mnsol.py) script using the [conda-lock_linux-64.yml](../../../data_generation/conda-lock_linux-64.yml) environment. 
Some ligands could not be charged with all methods, the following lists the ligands that could not be charged with each method:

- `am1bcc_oe`, `am1bcc_at`, `am1bccelf10_oe`: hydrogen
- `nagl_off`: hydrogen, tetramethylsilane

## Subsets

Regenerate all subsets: run `run.sh` in `data_generation/`.

### subset_openff_filtered.json
1587 systems retained from ~3037 MNSol database entries, covering 232 unique solutes across multiple solvent environments. Excluded: 219 charged systems; 528 entries with solvents on the skip list; 313 entries with solvents bearing disqualifying SMIRKS; 198 with missing solvent SMILES; 44 with solutes containing out-of-scope elements; 42 self-solvation entries (solute = solvent); 40 with solvents containing out-of-scope elements; 24 with solutes on the skip list; 18 with undefined solvent stereochemistry; 10 with disqualifying solute SMIRKS; 10 with undefined solute stereochemistry; and 4 with missing solute SMILES. Entries may satisfy multiple exclusion criteria; totals are not additive.

### subset_openff_small.json
400 systems, 188 unique solutes, with a per-solute cap of 3. Solute repetition: 67 appear once (59 pool-exhausted, 8 selection-limited), 30 appear twice (22 pool-exhausted, 8 selection-limited), and 91 appear three times (13 pool-exhausted, 78 at cap). Selection seeded the 189 MNSol–FreeSolv overlap solutes, then applied ChemicalEnvironment coverage fill and UMAP-guided diversification with 500 rebalancing swaps (483 same-solute, 17 cross-solute). No additional solutes were needed for coverage fill beyond the seeded overlap. 183 of 188 selected solutes (97.3%) overlap with the FreeSolv filtered pool. 25 ChemicalEnvironments are absent from the filtered pool and unrepresented.