1. Protein structure was prepared following the maestro instructions: https://github.com/OpenFreeEnergy/IndustryBenchmarks2024/pull/89.
2.The structure did not contain any cofactors, the ions were renamed from X1 -> Ca, and no non
natural amino acids were found.
2. Ligands were left unchanged, no protonation or conformational states were detected.
3. tested successfully the input protein with `input_validation.py`

using maestro v2024-2

Done by H. Baumann
* Added caps to all charged termini, except the N terminus in chain A which was left charge to form a salt bridge with an ASP to stabilize the binding site.
* Used the `charge_molecules.py` script to assign charges to the ligands
* Used the `generate_lomap_networks.py` script to generate the networks
* Use the `generate_ross_ref_data.py` script to generate the experimental reference data json file
