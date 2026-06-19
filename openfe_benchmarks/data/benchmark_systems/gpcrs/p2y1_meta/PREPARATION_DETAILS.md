- Prepared solvated protein-membrane complex was taken from https://github.com/schrodinger/public_binding_free_energy_benchmark/tree/main/fep_benchmark_inputs/fep_plus_inputs/gpcrs/p2y1_meta_sub.fmp
- The .fmp file was loaded into maestro and the lipids, solvent, and protein were saved separately as PDB files.
- The individual components were then combined into a single PDB file using openmm
- Box vectors were taken from the .fmp file and stored in the CRYST1 record in the PDB file
- Ligands were taken from https://github.com/schrodinger/public_binding_free_energy_benchmark/tree/main/fep_benchmark_inputs/structure_inputs/gpcrs/p2y1_meta_sub_ligands.sdf

Lomap network was generated using the generate_lomap_networks.py

Ligand partial charges were obtained using the charge_molecules.py script in the conda-lock environment.
