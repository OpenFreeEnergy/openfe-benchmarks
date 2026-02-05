# OpenFE Benchmarks Data

This directory contains benchmark system data for OpenFE (Open Free Energy) calculations, including protein structures, ligands, network, and remediated industry benchmark systems.

## Industry Benchmark Systems

### Structure

The `benchmark_systems/` directory contains remediated benchmark inputs organized as:

```
benchmark_systems/
├── <benchmark_set>/           # e.g., jacs_set, fragments, janssen_bace, mcs_docking_set
│   └── <system_name>/         # e.g., p38, tyk2, mcl1
│       ├── PREPARATION_DETAILS.md          # System preparation documentation generated
│       ├── protein.pdb        # Protein structure with cocrystallized waters & ions
│       ├── ligands.sdf        # Ligands without charges
│       ├── cofactors.sdf                   # (Optional) System cofactors without charges
│       ├── ligands_<charge_type>.sdf       # Ligands with specified partial charges generated with charge_molecules.py
│       ├── cofactors_<charge_type>.sdf     # (Optional) System cofactors with charges generated with charge_molecules.py
│       └── <network_name>.json             # Network mappings (e.g., industry_benchmarks_network.json) generated with generate_industry_lomap_networks.py
```

**Notes**:
- All charges generated using `openfe_benchmarks/data/data_generation/conda-lock_linux-64.yml`.

### Using Industry Benchmark Systems

The module provides a Python API for accessing benchmark systems:

```python
from openfe_benchmarks.data import (
    get_benchmark_data_system,
    list_benchmark_sets,
    list_data_systems,
    BenchmarkData,
    PARTIAL_CHARGE_TYPES,
    get_benchmark_set_data_systems,
)

# Discover available benchmark sets
benchmark_sets = list_benchmark_sets()
print(f"Available sets: {benchmark_sets}")

# List systems in a benchmark set
systems = list_data_systems('jacs_set')
print(f"Systems in jacs_set: {systems}")

# Load a specific benchmark system
system = get_benchmark_data_system('jacs_set', 'p38')

# Load all benchmark systems in a set
systems = get_benchmark_set_data_systems('jacs_set')

# Access system components
print(f"Protein: {system.protein}")
print(f"Ligands: {system.ligands}")  # Dict mapping charge type to file path
print(f"Cofactors: {system.cofactors}")  # Dict mapping charge type to file path (may be empty)
print(f"Networks: {system.network}")  # Network JSON file paths
```

For a complete tutorial on interaction with benchmark dataset, see `examples/1_initializing_benchmark_data.ipynb`.
For an example on applying a dataset, see `examples/2_benchmark_data_with_openfe.ipynb`.

## Additional Resources

- **OpenFE Documentation**: https://docs.openfree.energy
- **GUFE Documentation**: https://gufe.openfree.energy
- **Example Notebooks**: See `examples/` directory for:
  - `1_initializing_benchmark_data.ipynb` - How to use the benchmark systems API
  - `building_networks.ipynb` - Creating ligand networks

## Contributing

When adding new benchmark systems:

1. Follow the directory structure in `benchmark_systems/`
2. Document preparation details in `PREPARATION_DETAILS.md` within the system directory
3. Name files according to conventions:
   - `protein.pdb` for protein structures
   - `ligands.sdf` for ligands without charges
   - `cofactors.sdf` for cofactors without charges (if present)
   - `PREPARATION_DETAILS.md` (required) - system preparation documentation
4. Generate network mapping files in the `.json` format with ``openfe_benchmarks/data/data_generation/generate_industry_lomap_networks.py``
   - `*network.json` for network files (JSON format)
5. Add new charge types to `PARTIAL_CHARGE_TYPES` in `data/__init__.py` if necessary
6. Generate `ligands_<charge_type>.sdf` files with ``openfe_benchmarks/data/data_generation/charge_molecules.py``
6. (If present) Generate `cofactors_<charge_type>.sdf` files with ``openfe_benchmarks/data/data_generation/charge_molecules.py``