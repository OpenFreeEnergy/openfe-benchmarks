# OpenFE Benchmarks Data

This directory contains benchmark system data for OpenFE (Open Free Energy) calculations, including protein structures, ligands, ligand networks, and sometimes cofactors for remediated industry benchmark systems.

## Industry Benchmark Systems

### Structure

The `benchmark_systems/` directory contains remediated benchmark inputs organized as:

```
benchmark_systems/
├── <benchmark_set>/                     # e.g., jacs_set, fragments, janssen_bace, mcs_docking_set
│   └── <system_name>/                   # e.g., p38, tyk2, mcl1
│       ├── PREPARATION_DETAILS.md       # (Required) System preparation documentation generated
│       ├── ligands.sdf                  # (Required) Ligands without charges
│       ├── protein.pdb                  # (Optional) Protein structure with cocrystallized waters & ions. Include `protein` tag.
│       ├── <network_name>.json          # (Optional) Ligand network mappings (e.g., industry_benchmarks_network.json)
│       ├── cofactors.sdf                # (Optional) System cofactors without charges.  Include `cofactor` tag.
│       ├── ligands_<charge_type>.sdf    # (Required) Ligands with specified partial charges generated with charge_molecules.py
│       └── cofactors_<charge_type>.sdf  # (Optional) System cofactors with charges generated with charge_molecules.py
```

**Notes**:
- All charges generated using `openfe_benchmarks/data/data_generation/conda-lock_linux-64.yml`.

### Using Industry Benchmark Systems

For a complete tutorial on interaction with benchmark dataset, see `examples/1_initializing_benchmark_data.ipynb`.
For an example on applying a dataset, see `examples/2_benchmark_data_with_openfe.ipynb`.

## Additional Resources

- **OpenFE Documentation**: https://docs.openfree.energy
- **GUFE Documentation**: https://gufe.openfree.energy
- **Example Notebooks**: See `examples/` directory for:
  - `1_initializing_benchmark_data.ipynb` - How to use the benchmark systems API
  - `2_benchmark_data_with_openfe.ipynb`  - Using a BenchmarkData system with OpenFE tools
  - `building_networks.ipynb` - Creating ligand networks

## Contributing

When adding new benchmark systems:

1. Follow the directory structure in `benchmark_systems/` described in the [Structure](#structure) section.
2. Document preparation details in `PREPARATION_DETAILS.md` within the system directory
3. Name files according to conventions:
   - `ligands.sdf` for ligands without charges
   - `protein.pdb` for protein structures (if present)
   - `cofactors.sdf` for cofactors without charges (if present)
   - `PREPARATION_DETAILS.md` (required) - system preparation documentation
4. Generate ligand network mapping files, `*network*.json` (JSON format)
5. (If necessary) Add new charge types to `PARTIAL_CHARGE_TYPES` in `data/__init__.py`
6. Generate `ligands_<charge_type>.sdf` files with ``openfe_benchmarks/data/data_generation/charge_molecules.py``
6. (If present) Generate `cofactors_<charge_type>.sdf` files with ``openfe_benchmarks/data/data_generation/charge_molecules.py``