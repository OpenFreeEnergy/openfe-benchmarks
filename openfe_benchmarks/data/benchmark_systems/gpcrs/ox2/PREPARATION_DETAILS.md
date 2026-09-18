- Prepared solvated protein-membrane complex was taken from https://github.com/schrodinger/public_binding_free_energy_benchmark/tree/main/fep_benchmark_inputs/fep_plus_inputs/gpcrs/ox2_hip_custcore.fmp
- The .fmp file was loaded into maestro and the lipids, solvent, and protein were saved separately as PDB files.
- The individual components were then combined into a single PDB file using openmm
- Box vectors were taken from the .fmp file and stored in the CRYST1 record in the PDB file
- Ligands were taken from https://github.com/schrodinger/public_binding_free_energy_benchmark/tree/main/fep_benchmark_inputs/structure_inputs/gpcrs/ox2_hip_custcore_ligands.sdf
- Ligand set1: 13 was energy minimized using openmm for 500 steps ("openff-2.2.0.offxml", gasteiger partial charges) since two hydrogens were so close that partial charge generation failed for this conformation

Lomap network was generated using the generate_lomap_networks.py

Ligand partial charges were obtained using the charge_molecules.py script in the conda-lock environment.
