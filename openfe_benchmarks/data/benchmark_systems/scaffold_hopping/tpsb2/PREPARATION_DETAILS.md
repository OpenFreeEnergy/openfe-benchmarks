# System preparation details

By H. Baumann

* Prep done with PyMol 2.5.2
* Input structures taken from https://github.com/OpenFreeEnergy/IndustryBenchmarks2024/blob/a84180505e97de3ee0de55176fd1009240e889c2/industry_benchmarks/input_structures/original_structures/scaffold_hopping
* ACE and NME caps added
* TER removed from terminal NME caps
* Extracted the cofactor from the PDB and put it into an SDF
* Used the `charge_offmol.py` script to assign charges to the ligands and cofactor 
* Used the `generate_lomap_networks.py` script to generate the networks
* Use the `generate_ross_ref_data.py` script to generate the experimental reference data json file
