#!/usr/bin/env python
"""Orchestrate metadata preparation for alchemical archive submission.

Uses the Python API directly to invoke the three-step workflow:
1. Pull alchemical network from Alchemiscale
2. Generate computational results
3. Generate submission metadata
"""

import logging
import sys
from pathlib import Path

# Set up path for imports
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "scripts"))
from utils import setup_file_logger

# Import the core functions (not the CLI wrappers)
from _tmp_alchemiscale_gather import run_gather
from generate_results import run_generate_results
from prepare_metadata_submission import process_network

logger = logging.getLogger(__name__)


if __name__ == "__main__":
    setup_file_logger("log.txt", level=logging.INFO, print_console=True)
    logger.info("Starting metadata preparation workflow")

    # Step 1: Gather network from Alchemiscale
    logger.info("=" * 70)
    logger.info("Step: Pulling alchemical network")
    logger.info("=" * 70)
    run_gather(
        network_key="AlchemicalNetwork-2dd5d032b0228c7474eda50d8e064c2d-openff-test-openff_2_3_0_tyk2",
        allow_partial=False,
        output=".",
    )
    logger.info("✓ Pulling alchemical network completed successfully\n")

    # Step 2: Generate computational results
    logger.info("=" * 70)
    logger.info("Step: Generating computational_results.json")
    logger.info("=" * 70)
    run_generate_results(
        archive=(
            Path(
                "AlchemicalNetwork-2dd5d032b0228c7474eda50d8e064c2d-openff-test-openff_2_3_0_tyk2.json.bz2"
            ),
        ),
        system_group="jacs_set",
        system_name="tyk2",
        output_dir=Path("."),
    )
    logger.info("✓ Generating computational_results.json completed successfully\n")

    # Step 3: Generate submission metadata
    logger.info("=" * 70)
    logger.info("Step: Generating submission metadata")
    logger.info("=" * 70)
    process_network(
        input_files="AlchemicalNetwork*.json.bz2",
        output_dir=Path("."),
        submission_id="2026-06-22-tyk2-alchemicalarchive-test",
        tags="rbfe,benchmark,openfe",
        author=["Jennifer A. Clark"],
        license="CC-BY-4.0",
        system_group="jacs_set",
        system_name="tyk2",
        summary_suffix=(
            "This subset of the JACS set is meant to provide an example of an alchemical archive submission"
            " and provide an indication of the variability in results. For scripts to generate this network: "
            "github.com/openforcefield/alchemical-benchmark-resources/submissions/2026_03_17_openff-2.3.0_jacs_tyk2/alchemiscale_submission.ipynb"
        ),
    )
    logger.info("✓ Generating submission metadata completed successfully\n")

    logger.info("=" * 70)
    logger.info("✓ Workflow complete!")
    logger.info("=" * 70)
    logger.info("Check files in this directory:")
    logger.info("  - submission.yaml")
    logger.info("  - zenodo_description.md (do not include in submission)")
