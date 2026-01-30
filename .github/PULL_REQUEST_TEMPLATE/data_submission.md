# Pull Request Template: Benchmark Data Submission

## Checklist

Please ensure the following criteria are met before submitting the PR:

- [ ] **Data Format**: All files are in the expected format as outlined in `openfe_benchmarks/data/industry_benchmark_systems/README.md`.
- [ ] A `ligand.sdf` file is present
- [ ] A `protein.pdb` file is present for binding free energy calculations (bfe)
- [ ] **Documentation**: A `PREPARATION_DETAILS.md` file is included in the data directory, describing:
  - The source of the data
  - The method used to generate the data
  - Any assumptions or limitations
- [ ] The `charge_molecules.py` file was run using the `conda-lock_linux-64.yml` env on `ligands.sdf` for all partial charge types
- [ ] (If appropriate) The `charge_molecules.py` file was run using the `conda-lock_linux-64.yml` env on `cofactors.sdf` for all partial charge types
- [ ] Produce network with `generate_industry_lomap_networks.py`

Thank you for contributing to the OpenFE Benchmarks!