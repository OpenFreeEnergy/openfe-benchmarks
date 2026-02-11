# OpenFE Benchmarks

Curated benchmark inputs and result submission conventions for OpenFE-based free-energy calculations. The repository provides remediated system inputs (proteins, ligands, ligand networks), example scripts/notebooks, and a machine-readable results schema.

## Contents

- `openfe_benchmarks/data/` — remediated benchmark inputs (ligands, proteins, networks, preparation notes).
- `openfe_benchmarks/scripts/` — utility scripts and validation helpers.
- `openfe_benchmarks/results/` — canonical result-submission layout and schema (see `results/README.md`).
- `examples/` — runnable notebooks showing common workflows.

## Quick start

1. Install the package. No dependencies are needed to access data.

   ```bash
   pip install -e .
   ```

2. Browse benchmarking systems with [data](openfe_benchmarks/data/README.md)
3. Set-up calculations with scripts.
4. Submit benchmarking [results](openfe_benchmarks/results/README.md)

## License

This project is MIT-licensed. See `LICENSE` for details.
