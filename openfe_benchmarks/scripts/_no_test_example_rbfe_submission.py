#!/usr/bin/env python
"""Orchestrate metadata preparation for alchemical archive submission.

Uses the Python API directly to invoke the three-step workflow:
1. Pull alchemical network from Alchemiscale
2. Generate computational results
3. Generate submission metadata
"""

import logging
from pathlib import Path

from openfe_benchmarks.scripts.utils import setup_file_logger

# Import the core functions (not the CLI wrappers)
from openfe_benchmarks.scripts._tmp_alchemiscale_gather import run_gather
from openfe_benchmarks.scripts.generate_results_archives import run_generate_results
from openfe_benchmarks.scripts.prepare_metadata_submission import process_network

logger = logging.getLogger(__name__)

NETWORK_KEY = "your_network_key"
SYSTEM_GROUP = "jacs_set"
SYSTEM_SET = "tyk2"
OUTPUT_DIR = "output"

SUBMISSION_ID = "2026-06-22-tyk2-alchemicalarchive-test"
DATE = "2026-06-22"
AUTHOR = "Jane Doe"
SUMMARY_SUFFIX = (
    "This subset of the JACS set is meant to provide an example of an alchemical archive submission"
    " and provide an indication of the variability in results. For scripts to generate this network: "
    "github.com/openforcefield/alchemical-benchmark-resources/submissions/2026_03_17_openff-2.3.0_jacs_tyk2/alchemiscale_submission.ipynb"
)
TAGS = "rbfe,benchmark,openfe"
OPENFE_VER = "1.8.0"
OPENMM_VER = "8.2.0"
OFFTOOL_VER = "0.18"

if __name__ == "__main__":
    setup_file_logger("log.txt", level=logging.INFO, print_console=True)
    logger.info("Starting metadata preparation workflow")

    # Step 1: Gather network from Alchemiscale
    logger.info("=" * 70)
    logger.info("Step: Pulling alchemical network")
    logger.info("=" * 70)
    run_gather(
        network_key=NETWORK_KEY,
        allow_partial=False,
        output=OUTPUT_DIR,
    )
    logger.info("✓ Pulling alchemical network completed successfully\n")

    # Step 2: Generate computational results
    logger.info("=" * 70)
    logger.info("Step: Generating computational_results.json.bz2")
    logger.info("=" * 70)
    run_generate_results(
        archive=Path(f"{OUTPUT_DIR}/{NETWORK_KEY}.json.bz2"),
        system_group=SYSTEM_GROUP,
        system_name=SYSTEM_SET,
        output_dir=Path(OUTPUT_DIR),
    )
    logger.info("✓ Generating computational_results.json.bz2 completed successfully\n")

    # Step 3: Generate submission metadata
    logger.info("=" * 70)
    logger.info("Step: Generating submission metadata")
    logger.info("=" * 70)
    process_network(
        input_files=f"{OUTPUT_DIR}/AlchemicalNetwork*.json.bz2",
        output_dir=Path(OUTPUT_DIR),
        submission_id=SUBMISSION_ID,
        tags=TAGS,
        author=[AUTHOR],
        license="CC-BY-4.0",
        system_group=SYSTEM_GROUP,
        system_name=SYSTEM_SET,
        submission_date=DATE,
        summary_suffix=SUMMARY_SUFFIX,
        openfe_version=OPENFE_VER,
        openmm_version=OPENMM_VER,
        openff_toolkit_version=OFFTOOL_VER,
    )
    logger.info("✓ Generating submission metadata completed successfully\n")

    logger.info("=" * 70)
    logger.info("✓ Workflow complete!")
    logger.info("=" * 70)
    logger.info("Check files in this directory:")
    logger.info("  - submission.yaml")
    logger.info("  - zenodo_description.md (do not include in submission)")
