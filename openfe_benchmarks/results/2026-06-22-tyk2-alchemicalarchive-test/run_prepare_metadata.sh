#!/bin/bash
# Generate submission metadata for all networks (charge_changes and jacs_set)
# This script processes all networks from the ResultSubmission folder

set -e  # Exit on error

echo "Pulling alchemical network (will be hosted on zenodo instead of submitted in this repository)"
echo ""

micromamba run -n openfe-benchmarks-test python \
    ../../scripts/_tmp_alchemiscale_gather.py \
    --network-key AlchemicalNetwork-2dd5d032b0228c7474eda50d8e064c2d-openff-test-openff_2_3_0_tyk2 \
    --no-allow-partial \
    --output "."

echo "Generating computational_results.json..."
echo ""

micromamba run -n openfe-benchmarks-test python \
    ../../scripts/example_generate_results.py \
    --archive "AlchemicalNetwork-2dd5d032b0228c7474eda50d8e064c2d-openff-test-openff_2_3_0_tyk2.json.bz2" \
    --system-group "jacs_set" \
    --system-name "tyk2" \
    --output_dir "."

echo "Generating submission metadata for all benchmark networks..."
echo ""

micromamba run -n openfe-benchmarks-test python \
    ../../scripts/prepare_metadata_submission.py \
    "AlchemicalNetwork*.json.bz2" \
    --system-group "jacs_set" \
    --system-name "tyk2" \
    --output-dir . \
    --submission-id "2026-06-22-tyk2-alchemicalarchive-test" \
    --tags "rbfe, benchmark, openfe" \
    --author "Jennifer A. Clark" \
    --license "CC-BY-4.0" \
    --summary-suffix "This subset of the JACS set is meant to provide an example of an alchemical archive submission and provide an indication of the variability in results."

echo ""
echo "Done! Check the files in this directory:"
echo "  - submission.yaml"
echo "  - zenodo_description.md"
