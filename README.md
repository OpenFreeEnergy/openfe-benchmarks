# openfe-benchmarks

Benchmark datasets and archived computational results for OpenFE benchmarking.

## Quick links

- [How to install](#how-to-install)
- [How to access benchmark results](#how-to-access-benchmark-results)
- [How to add a benchmark result](#how-to-add-a-benchmark-result)
- [How to add a dataset](#how-to-add-a-dataset)

## How to install

Install in a Python 3.11+ environment:

```bash
git clone https://github.com/OpenFreeEnergy/openfe-benchmarks.git
cd openfe-benchmarks
pip install -e .
```

If you use the local OpenFE conda setup, run commands with:

```bash
micromamba run -n openfe <command>
```

## How to access benchmark results

Use the results API to load one submission by ID, or filter many submissions.

```python
from openfe_benchmarks.results import get_benchmark_results, filter_results

result = get_benchmark_results("2026-03-18-openmm-840-qa-testing") # from submission id
print(result.title)
print(result.raw_results.keys())  # dg / ddg

rbfe_submissions = filter_results(calculation_type="rbfe", tags="jacs_set")
print(len(rbfe_submissions))
```

Submission layout is stored in:

- openfe_benchmarks/results/<submission_id>/submission.yaml
- openfe_benchmarks/results/<submission_id>/computational_results.json.bz2

Examples for analysis and plotting:

- examples/4_benchmark_result_plot.ipynb
- examples/5_multicomparison_plots.ipynb
- openfe_benchmarks/scripts/_example_plot_rbfe.py
- openfe_benchmarks/scripts/_no_test_example_plot_asfe.py

## How to add a benchmark result

Goal: add one new submission directory under openfe_benchmarks/results/.

Use this 3-step flow:

1. Gather archive(s) from Alchemiscale or local output.
2. Build computational results.
3. Build submission metadata.

Repository examples that orchestrate the full metadata workflow:

- openfe_benchmarks/scripts/_no_test_example_rbfe_asfe_submission.py
- openfe_benchmarks/scripts/_no_test_example_mutlinetwork_rbfe_submission.py

Core scripts:

- openfe_benchmarks/scripts/generate_results_archives.py
- openfe_benchmarks/scripts/prepare_metadata_submission.py

Minimal example commands:

```bash
# 1) Create computational_results.json.bz2 from one archive
python openfe_benchmarks/scripts/generate_results_archives.py \
  --archive path/to/AlchemicalNetwork-<hash>.json.bz2 \
  --output-dir openfe_benchmarks/results/<submission_id> \
  --system-group jacs_set \
  --system-name tyk2

# 2) Create submission.yaml and zenodo_description.md
python openfe_benchmarks/scripts/prepare_metadata_submission.py \
  path/to/AlchemicalNetwork-<hash>.json.bz2 \
  --output-dir openfe_benchmarks/results/<submission_id> \
  --submission-id <submission_id> \
  --submission-date YYYY-MM-DD \
  --author "Your Name" \
  --results-file computational_results.json.bz2
```

## How to add a dataset

Goal: add a new benchmark system under openfe_benchmarks/data/benchmark_systems/.

1. Create the dataset folder at:
   openfe_benchmarks/data/benchmark_systems/<system_group>/<system_name>/
2. Add required inputs:
   - PREPARATION_DETAILS.md
   - ligands.sdf
   - ligands_<charge_type>.sdf
3. Add optional files as needed:
   - protein.pdb
   - cofactors.sdf
   - cofactors_<charge_type>.sdf
   - network JSON files
   - experimental data JSON (for bfe/sfe benchmarks)
4. Update benchmark index tags in:
   openfe_benchmarks/data/benchmark_system_indexing.yml
5. Validate data access with notebooks and API examples.

Helpful references:

- openfe_benchmarks/data/README.md
- examples/1_initializing_benchmark_data.ipynb
- examples/2_benchmark_data_with_openfe.ipynb
- examples/building_networks.ipynb

## External docs

- OpenFE docs: https://docs.openfree.energy
- GUFE docs: https://gufe.openfree.energy