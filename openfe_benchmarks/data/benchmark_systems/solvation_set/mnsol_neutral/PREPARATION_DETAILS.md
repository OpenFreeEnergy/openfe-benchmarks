Done by Jennifer Clark on 2026-02-25

Systems extracted from the MNSol database (https://doi.org/10.13020/3eks-j059).

## Filtering

The reference data was generated using the [generate_mnsol_data.py](../../../data_generation/generate_mnsol_data.py) script. Entries were excluded if:

- The solute or solvent name was not present in `mnsol-name-to-smiles.json`
- Molecules with ambiguous isomeric structure were excluded including: 'bromotoluene', 'chlorotoluene', 'dichloroethane', 'fluoroctane', 'trimethylbenzene'
- The solute is a radical (name contains "radical")
- The solute is a dimer (MNSol Level1=14)
- The net charge is non-zero

Conformer generation failed for: 6,7,8,9,10,10-hexachloro-1,5,5a,6,9,9a-hexahydro-6,9-methano-2,4,3-benzodioxathiepine-3-oxide(endosulfanalpha)

After filtering, **2468 systems** are retained.

## Charging Solutes / Solvents

Charges were generated using the [charge_mnsol.py](../../../data_generation/charge_mnsol.py) script using the [conda-lock_linux-64.yml](../../../data_generation/conda-lock_linux-64.yml) environment. 
Some ligands could not be charged with all methods, the following lists the ligands that could not be charged with each method:

`am1bcc_oe`: "hydrogen"
`am1bcc_at`: "hydrogen"
`am1bccelf10_oe`: "hydrogen"
`nagl_off`: "hydrogen", "tetramethylsilane"

## Subsets

Of the 3037 systems imported from MNSol, two subsets were made:

- "subset_openff_filtered": 1587 
- "subset_openff_small": 581 systems were down selected from "subset_openff_filtered" using Morgan fingerprint Tanimoto distance and specifying inclusion of checkmol functional groups.

Rows are excluded if the solute or solvent is water, absent from ``mnsol-name-to-smiles.json``, contains 'radical', is an explicit-solvent entry (Level1==14), is charged, is a self-solvation pair, matches a disqualifying SMIRKS (long chains, 1,3-dicarbonyls, or dissociating molecules, i.e., HBr or HCl), has undefined stereochemistry, or contains elements outside ``ALLOWED_ELEMENTS``.

## Notes

- Licensing prevented including MNSol experimental values. After properly obtaining the dataset, a user can locally generate the needed experimental data file with the [generate_mnsol_data.py](../../../data_generation/generate_mnsol_data.py) script.
- Experimental uncertainties are set to 0.2 kcal/mol for all neutral entries, following the recommendation in the MNSol documentation.