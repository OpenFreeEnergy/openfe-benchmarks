# Open Free Energy Benchmark Systems

A set of benchmark systems to validate the OpenFE components.

## Contents

The `openfe_benchmark` repository contains:
  * A file for each system in the benchmark set (e.g. `tyk2.py`). These hold
    methods for creating the free energy network and the system components
    necessary to calculate the benchmark.
  * A `data` folder with a set of PDB files for each host system and SDF files
    for each set of ligands for each system.
  * `util.py` a set of utility methods and base classes used throughout this
    repository.
  * `examples` a directory with notebooks for example calculations and analyses

## Relative binding free energies

The following systems are obtained from the [OpenFF Protein-Ligand Benchmarks](https://github.com/openforcefield/protein-ligand-benchmark):
  * `TYK2`
  * `PTP1B`
