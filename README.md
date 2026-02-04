# Open Free Energy Benchmark Systems

A set of benchmark systems to validate the OpenFE components.

## Contents

The `openfe_benchmark` repository contains:
  * A file for each system in the benchmark set (e.g. `tyk2.py`). These hold
    methods for creating the free energy network and the system components
    necessary to calculate the benchmark.
  * A [`data`](/Users/jenniferclark/bin/openfe-benchmarks/openfe_benchmarks/data/README.md) directory with a set of PDB files for each host system and SDF files
    for each set of ligands for each system.
  * A `scripts` directory with examples of how to set up calculations such as:
   - bfe: binding free energy
   - hfe: hydration free energy
   - sfe: solvation free energy (requires external tools)
  * `examples` a directory with notebooks for example calculations and analyses

## Available Systems

The available systems are indexed in the [benchmark_system_indexing.yml](openfe_benchmarks/data/benchmark_system_indexing.yml) file. This file provides a comprehensive list of the benchmark systems along with tags representing whether they can be used for `bfe` or `sfe` calculations and whether there are cofactors in the system. These systems can also be explored in Python, as demonstrated in the notebook [`1_initializing_benchmark_data.ipynb`](examples/1_initializing_benchmark_data.ipynb).

## Additional Resources

- **OpenFE Documentation**: https://docs.openfree.energy
- **GUFE Documentation**: https://gufe.openfree.energy
- **Example Notebooks**: See `examples/` directory for:
  - `using_benchmark_systems.ipynb` - How to use the benchmark systems API
  - `building_networks.ipynb` - Creating ligand networks