# Benchmark Results Submission

## Description
<!-- Systems calculated, single experiment, notable settings. -->


## Change level:
See [RELEASE.md](../RELEASE.md).
- [ ] `major`: changes existing data
- [ ] `minor`: additive
- [ ] `bug`: fix, script, docs, new results

## Submission
- [ ] Labeled PR with version level `major`, `minor`, or `bug`
- [ ] `submission.yaml` made with `openfe_benchmarks/scripts/prepare_metadata_submission.py`
- [ ] `computational_results.json.bz2` made with `openfe_benchmarks/scripts/generate_results_archives.py`
- [ ] Single experiment; settings differing from OpenFE defaults described in `summary` of this PR submission
- [ ] [CHANGELOG.md](../../CHANGELOG.md) updated under `Unreleased`

## Checked by CI
- `submission.yaml` is valid and loads
- `submission_id` matches directory name and is unique
- Results file exists at the path in `results`

## Notes
