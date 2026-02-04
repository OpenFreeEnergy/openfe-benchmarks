# Pull Request Template: Benchmark Data Submission

Only one of the following check lists need be used.

## New System Checklist

Please ensure the following criteria are met before submitting the PR:

- [ ] **Data Format**: All files are in the expected format as outlined in `openfe_benchmarks/data/benchmark_systems/README.md`.
- [ ] A `ligand.sdf` file is present
- [ ] The `charge_molecules.py` file was run using the `conda-lock_linux-64.yml` env on `ligands.sdf`
- [ ] (If appropriate) A `protein.pdb` file is present for binding free energy calculations (bfe)
- [ ] (If appropriate) The `charge_molecules.py` file was run using the `conda-lock_linux-64.yml` env on `cofactors.sdf`
- [ ] (If appropriate) Produce ligand network for r*fe calculations
- [ ] **Documentation**: A `PREPARATION_DETAILS.md` file is included in the data directory, describing:
  - The source of the data
  - The method used to generate the data
  - Any assumptions or limitations
- [ ] Ensure this new set/system is represented in `benchmark_system_indexing.yml`
- [ ] If there is a new partial charge method, it must be added to PARTIAL_CHARGE_TYPES in `_benchmark_systems.py`

## Update System Checklist

Please ensure the following criteria are met before submitting the PR:

- [ ] **Data Format**: All files are in the expected format as outlined in `openfe_benchmarks/data/benchmark_systems/README.md`.
- [ ] (If appropriate) The `charge_molecules.py` file was run using the `conda-lock_linux-64.yml` env on `ligands.sdf`
- [ ] (If appropriate) The `charge_molecules.py` file was run using the `conda-lock_linux-64.yml` env on `cofactors.sdf`
- [ ] (If appropriate) Produce ligand network for r*fe calculations
- [ ] **Documentation**: A `PREPARATION_DETAILS.md` file includes a changelog of date and what was changed in the system
- [ ] If there is a new partial charge method, it must be added to PARTIAL_CHARGE_TYPES in `_benchmark_systems.py`

Thank you for contributing to the OpenFE Benchmarks!