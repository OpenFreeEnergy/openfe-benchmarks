# Pull Request Template: Benchmark Results Submission

## Description
[Provide a brief description of this benchmark submission, including what systems were calculated and any notable features]

## Required Files Checklist

Please ensure the following files are present in your submission directory:

- [ ] **`submission.yaml`**: Complete metadata file with all required fields using `openfe_benchmarks/scripts/prepare_metadata_submission.py`
- [ ] **`computational_results.json.bz2`**: Compressed results file using `openfe_benchmarks/scripts/generate_results_archives.py`
- [ ] **Archive DOI**: Results uploaded to long-term archive (Zenodo, etc.) with DOI included in submission.yaml

## Validation Checklist

- [ ] **YAML validation**: `submission.yaml` is valid YAML and loads without errors
- [ ] **ID consistency**: `submission_id` in YAML matches directory name
- [ ] **Results file exists**: Compressed results file exists at path specified in `results` field
- [ ] **Archive accessible**: DOI resolves and archive is publicly accessible
- [ ] **Network keys valid**: AlchemicalNetwork keys in `benchmark_data` are valid
- [ ] **Calculation type valid**: Type matches actual calculations performed
- [ ] **No duplicate submission_id**: This submission_id is unique in the repository

## Testing

- [ ] CI validation checks pass

## Additional Notes

[Any additional context, special considerations, or notes about this submission]

---

Thank you for contributing to the OpenFE Benchmarks!