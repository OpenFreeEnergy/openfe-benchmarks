#!/usr/bin/env python
"""Orchestrate metadata preparation for alchemical archive submission.

Uses the Python API directly to invoke the three-step workflow:
1. Pull alchemical network from Alchemiscale
2. Generate computational results
3. Generate submission metadata
"""

import os
import logging
from pathlib import Path

from openfe_benchmarks.scripts.utils import setup_file_logger

# Import the core functions (not the CLI wrappers)
from openfe_benchmarks.scripts._tmp_alchemiscale_gather import run_gather
from openfe_benchmarks.scripts.generate_results_archives import run_generate_results
from openfe_benchmarks.scripts.prepare_metadata_submission import process_network

logger = logging.getLogger(__name__)

SYSTEM_GROUP_SET_KEY = (
    [
        "jacs_set",
        "bace",
        "AlchemicalNetwork-58ef1baba714145eaddb8f81abf609ce-openff-openff_3_0_0_alpha1b_opc3-rbfe_pontibus_jacs_set_bace",
    ],
    [
        "jacs_set",
        "cdk2",
        "AlchemicalNetwork-8b3e6b25fe046924998a1f520c798b19-openff-openff_3_0_0_alpha1b_opc3-rbfe_pontibus_jacs_set_cdk2",
    ],
    [
        "jacs_set",
        "jnk1",
        "AlchemicalNetwork-6cea205bc124349b8add8824e2060147-openff-openff_3_0_0_alpha1b_opc3-rbfe_pontibus_jacs_set_jnk1",
    ],
    [
        "jacs_set",
        "mcl1",
        "AlchemicalNetwork-ee5422ba932ced782b51da9bb9eb24ff-openff-openff_3_0_0_alpha1b_opc3-rbfe_pontibus_jacs_set_mcl1",
    ],
    [
        "jacs_set",
        "p38",
        "AlchemicalNetwork-c613a533bdaa47f92e9a1cf4d58fe6dd-openff-openff_3_0_0_alpha1b_opc3-rbfe_pontibus_jacs_set_p38",
    ],
    [
        "jacs_set",
        "ptp1b",
        "AlchemicalNetwork-3885e858e4b36849e110edbebe348ae2-openff-openff_3_0_0_alpha1b_opc3-rbfe_pontibus_jacs_set_ptp1b",
    ],
    [
        "jacs_set",
        "thrombin",
        "AlchemicalNetwork-7fcf1d4cb0a6af702747de83fc23af69-openff-openff_3_0_0_alpha1b_opc3-rbfe_pontibus_jacs_set_thrombin",
    ],
    [
        "jacs_set",
        "tyk2",
        "AlchemicalNetwork-51caad3a5d3daea240406eb23fa6dad6-openff-openff_3_0_0_alpha1b_opc3-rbfe_pontibus_jacs_set_tyk2",
    ],
)
OUTPUT_DIR = "output"

SUBMISSION_ID = "2026-08-04-openff3.0.0-alpha1b_opc3-jacs"
DATE = "2026-08-04"
AUTHOR = "Jennifer A. Clark"
SUMMARY_SUFFIX = (
    "For scripts to generate this network: "
    "github.com/openforcefield/alchemical-benchmark-resources/submissions/2026_07_15_openff-3.0.0-alpha1b_opc3/alchemiscale_submission.ipynb"
)
TAGS = "rbfe,benchmark,openfe"
SMALL_MOL_FF = "openff-3.0.0-alpha1b"
WATER_MODEL = "opc3.offxml"
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
    for _, _, network_key in SYSTEM_GROUP_SET_KEY:
        if not os.path.isfile(os.path.join(OUTPUT_DIR, f"{network_key}.json.bz2")):
            run_gather(
                network_key=network_key,
                allow_partial=False,
                output=OUTPUT_DIR,
            )
    logger.info("✓ Pulling alchemical network completed successfully\n")

    # Step 2: Generate computational results
    logger.info("=" * 70)
    logger.info("Step: Generating computational_results.json")
    logger.info("=" * 70)
    systems_local = [
        [group, sys_set, os.path.join(OUTPUT_DIR, f"{network_key}.json.bz2")]
        for group, sys_set, network_key in SYSTEM_GROUP_SET_KEY
    ]
    run_generate_results(
        systems=systems_local,
        output_dir=Path(OUTPUT_DIR),
    )
    logger.info("✓ Generating computational_results.json completed successfully\n")

    # Step 3: Generate submission metadata
    logger.info("=" * 70)
    logger.info("Step: Generating submission metadata")
    logger.info("=" * 70)
    process_network(
        systems=systems_local,
        output_dir=Path(OUTPUT_DIR),
        submission_id=SUBMISSION_ID,
        tags=TAGS,
        author=[AUTHOR],
        license="CC-BY-4.0",
        submission_date=DATE,
        summary_suffix=SUMMARY_SUFFIX,
        small_molecule_forcefield=SMALL_MOL_FF,
        forcefields=[SMALL_MOL_FF, WATER_MODEL],
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
    logger.info("  - zenodo_description.txt (do not include in submission)")
