# Release Policy of OpenFE-Benchmarks

OpenFE-Benchmarks uses **Data-SemVer**: semantic versioning applied to published benchmark data, where "breaking" means breaking for someone consuming the data.

- Version format is `X.Y.Z`, which is PEP 440 compliant.
- Releases are git tags `vX.Y.Z` on `main`. The package version is derived from the tag by `setuptools-scm`.

## Version levels

| Level | Meaning | Examples |
|---|---|---|
| **X** (major) | Changes to existing data; not backward compatible | Editing or removing values in an existing benchmark data system; renaming or removing loader API or `submission.yaml` fields |
| **Y** (minor) | Independent additions; backward compatible | New datasets or systems; new charge sets; new results submission; new optional fields or script features |
| **Z** (bug) | Fixes and changes that do not alter published data | Script fixes, docs, CIs |

Tie-breaker: if a change fits more than one level, use the highest.

## Process

1. The PR author picks the level and labels the PR `major`, `minor`, or `bug`.
2. The PR author adds an entry under `Unreleased` in [CHANGELOG.md](CHANGELOG.md).
3. On a quarterly basis maintainers cut a release by moving `Unreleased` into a dated `X.Y.Z` section and tagging `vX.Y.Z` on `main`.

Note: Major changes may be submitted as PRs but are not merged without agreement between the Science Team Leads for both OpenFE and OpenFF.

## Release Cadence
Releases are made quarterly.

## Zenodo and citation

Citation and DOI are obtained from Zenodo / GitHub integration for each `vX.Y.Z` tag.
