# System preparation details

By J.T.Horton

* Prep done with PyMol 3.0.0
* Input structures taken from https://github.com/OpenFreeEnergy/IndustryBenchmarks2024/blob/a84180505e97de3ee0de55176fd1009240e889c2/industry_benchmarks/input_structures/original_structures/opls_stress
* ACE and NME caps added to chain A
* TER removed from NME cap
* Used the `charge_molecules.py` script to assign charges to the ligands 
* Used the `generate_lomap_networks.py` script to generate the networks
* Use the `generate_ross_ref_data.py` script to generate the experimental reference data json file