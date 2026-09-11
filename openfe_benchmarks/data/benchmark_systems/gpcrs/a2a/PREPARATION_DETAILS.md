- Prepared solvated protein-membrane complex was taken from https://github.com/schrodinger/public_binding_free_energy_benchmark/tree/main/fep_benchmark_inputs/fep_plus_inputs/gpcrs/a2a_hip278.fmp
- The .fmp file was loaded into maestro and the lipids, solvent, and protein were saved separately as PDB files.
- The individual components were then combined into a single PDB file using openmm
- Box vectors were taken from the .fmp file and stored in the CRYST1 record in the PDB file
- Ligands were taken from https://github.com/schrodinger/public_binding_free_energy_benchmark/tree/main/fep_benchmark_inputs/structure_inputs/gpcrs/a2a_hip278_ligands.sdf
- For ligands with multiple binding modes, the binding mode was retained that Ross et al. found to be more potent (present in this table https://github.com/schrodinger/public_binding_free_energy_benchmark/blob/main/21_4_results/ligand_predictions/gpcrs/a2a_hip278_sbpkacorr_out.csv)
- For ligands with multiple protonation states, the neutral form was retained

Lomap network was generated using the generate_lomap_networks.py

Ligand partial charges were obtained using the charge_molecules.py script in the conda-lock environment.
