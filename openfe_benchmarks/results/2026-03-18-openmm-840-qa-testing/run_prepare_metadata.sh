#!/bin/bash
# Generate submission metadata for all networks (charge_changes and jacs_set)
# This script processes all networks from the ResultSubmission folder

set -e  # Exit on error

echo "Generating submission metadata for all benchmark networks..."
echo ""

micromamba run -n openfe-benchmarks python \
    /Users/jenniferclark/bin/openfe-benchmarks/openfe_benchmarks/scripts/prepare_metadata_submission.py \
    "/Users/jenniferclark/OMSF/OpenFE/BenchmarkRepo/ResultSubmission/networks/*/*/*alchemicalnetwork.json" \
    --output-dir . \
    --submission-id "2026-03-18-openmm-840-qa-testing" \
    --tags "charge_change, rbfe, benchmark, openfe, openmm-840" \
    --author "Josh Horton" \
    --license "CC-BY-4.0" \
    --no-alchemiscale \
    --summary-suffix "Note this means the charge annihilation sets are not complete compared to what is in that system and should not be compared to other complete runs due to the missing edges."

echo ""
echo "Done! Check the files in this directory:"
echo "  - submission.yaml"
echo "  - zenodo_description.md"
