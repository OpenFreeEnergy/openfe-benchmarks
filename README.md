# Open Free Energy Benchmark Systems

A set of benchmark systems to validate the OpenFE components.

## Contents

Each directory contains:
  * `protein.pdb` a PDB of the biological unit components of the system
  * `ligands.sdf` an SDF containing all the small molecule components which
    undergo transformations
  * `system.py` contains `get_components` a method which returns a ligand
    transformation network, and both the protein and solvent components for
    the transformations of interest

## Relative binding free energies

The following systems are obtained from the [OpenFF Protein-Ligand Benchmarks](https://github.com/openforcefield/protein-ligand-benchmark):
  * `tyk2`
