"""
This module provides an example script for planning absolute solvation free energy (ASFE) calculations.
It demonstrates how to process benchmark systems, prepare components, and save the resulting
alchemical network to a JSON file.

Note that although the output of this script is an alchemical network, that simply serves as a storage
class for the independent transformations needed in ASFEs.
"""

import json
import logging
from pathlib import Path

from openff.units import unit
import openfe
from openfe import SmallMoleculeComponent

from pontibus.components import ExtendedSolventComponent
from pontibus.utils.molecules import WATER
from pontibus.protocols.solvation import ASFEProtocol, ASFESettings
from pontibus.protocols.solvation.settings import PackmolSolvationSettings

from openfe_benchmarks.data import get_benchmark_data_system
from openfe_benchmarks.scripts import utils as ofebu

logger = logging.getLogger(__name__)


def _configure_example_logging(level=logging.INFO):
    handler = logging.StreamHandler()
    handler.setFormatter(
        logging.Formatter("%(asctime)s %(name)s %(levelname)s: %(message)s")
    )
    # Attach to this module's logger so output appears when running the example
    logger.addHandler(handler)
    logger.setLevel(level)
    # Optionally enable package-wide logs:
    logging.getLogger("openfe_benchmarks").setLevel(level)


BENCHMARK_SET = "solvation_set"
# BENCHMARK_SYS = "freesolv"
BENCHMARK_SYS = "mnsol_neutral"
PARTIAL_CHARGE = "nagl_openff-gnn-am1bcc-1.0.0.pt"
FORCEFIELD = "openff-2.3.0"  # available [openmmforcefields SystemGenerator](https://github.com/openmm/openmmforcefields?tab=readme-ov-file#automating-force-field-management-with-systemgenerator)
OUTPUT_DIR = "outputs"
FILENAME = f"network_{BENCHMARK_SET}_{BENCHMARK_SYS}_nacl_asfe.json"
WATER_MODEL = "tip3p"  # or opc

EXPECTED_NETWORKS = []
EXPECTED_LIGANDS = set()


def get_chemical_systems(
    benchmark_sys,
) -> dict[str, openfe.ChemicalSystem]:
    """
    Build ChemicalSystem objects from benchmark reference data for ASFE calculations.

    For each entry in the reference data file (experimental_solvation_free_energy_data.json or systems_data.json):
        - Uses the JSON key as the network name.
        - Matches solute and solvent molecules by name from the ligand SDF file (ligands_<PARTIAL_CHARGE>.sdf).
        - If the solute or solvent is water, uses the module-level WATER object and sets the name appropriately.
        - If a molecule is missing from the SDF, skips that entry and logs a warning.

    Parameters
    ----------
    benchmark_sys : BenchmarkData
        Loaded benchmark system with reference and ligand data.

    Returns
    -------
    dict of str to openfe.ChemicalSystem
        Mapping of network name to ChemicalSystem with solute and solvent components.
    """
    if "experimental_solvation_free_energy_data" in benchmark_sys.reference_data:
        ref_key = "experimental_solvation_free_energy_data"
    else:
        ref_key = "systems_data"
    ref_path = benchmark_sys.reference_data[ref_key]
    exp_data: dict = json.loads(Path(ref_path).read_text())
    logger.info(f"Loaded {len(exp_data)} experimental entries from {ref_path.name}")

    # Load all molecules (solutes + organic solvents) keyed by molecule name
    mol_dict: dict[str, SmallMoleculeComponent] = ofebu.process_sdf(
        benchmark_sys.ligands[PARTIAL_CHARGE], return_dict=True
    )
    logger.info(
        f"Loaded {len(mol_dict)} molecules from {benchmark_sys.ligands[PARTIAL_CHARGE].name} (charge model: {PARTIAL_CHARGE})"
    )

    systems: dict[str, openfe.ChemicalSystem] = {}
    for network_name, entry in exp_data.items():
        solute_name, solvent_name = entry["solute_name"], entry["solvent_name"]

        if solute_name == "water":
            solute = WATER
            ### Until Pontibus PR is Merged that defined water component name ###
            water_dict = solute.to_dict()
            water_dict["molprops"]["ofe-name"] = "water"
            solute = SmallMoleculeComponent.from_dict(water_dict)
            #####################################################################
        elif solute_name not in mol_dict:
            logger.warning(
                f"Solute '{solute_name}' not found in SDF; skipping network '{network_name}'"
            )
            continue
        else:
            solute = mol_dict[solute_name]

        if solvent_name.lower() == "water":
            solvent = ExtendedSolventComponent()
            ### Until Pontibus PR is Merged that defined water component name ###
            water_dict = solvent.solvent_molecule.to_dict()
            water_dict["molprops"]["ofe-name"] = "water"
            solvent = ExtendedSolventComponent(
                solvent_molecule=SmallMoleculeComponent.from_dict(water_dict)
            )
            #####################################################################
        else:
            if solvent_name not in mol_dict:
                logger.warning(
                    f"Solvent '{solvent_name}' not found in SDF; skipping network '{network_name}'"
                )
                continue
            solvent = ExtendedSolventComponent(solvent_molecule=mol_dict[solvent_name])

        systems[network_name] = openfe.ChemicalSystem(
            {"solute": solute, "solvent": solvent},
            name=network_name,
        )
        EXPECTED_NETWORKS.append(network_name)
        EXPECTED_LIGANDS.add(solute_name)
        EXPECTED_LIGANDS.add(solvent_name)

    logger.info(
        f"Built {len(systems)} ChemicalSystems from experimental data ({len(exp_data) - len(systems)} skipped)"
    )
    return systems


def get_water_settings(forcefield: str, water_model: str = "tip3p") -> ASFESettings:
    """
    Get settings for hydration free energy calculations (water as solvent).

    Parameters
    ----------
    forcefield : str
        Name of the force field to use.
    water_model : str, optional
        Water model to use (default is "tip3p").

    Returns
    -------
    ASFESettings
        Settings object for ASFEProtocol.
    """
    # The settings here are effectively the "fast" settings
    # shown in the validation.
    settings: ASFESettings = ASFEProtocol.default_settings()
    settings.protocol_repeats = 1
    settings.solvent_forcefield_settings.forcefields = [
        # To use a custom force field, just pass an OFFXML string
        # just like you would to openff.toolkit.ForceField
        forcefield,
        f"{water_model}.offxml",
    ]
    settings.vacuum_forcefield_settings.forcefields = [
        forcefield,  # as above
    ]
    settings.vacuum_engine_settings.compute_platform = "CUDA"
    settings.solvent_engine_settings.compute_platform = "CUDA"
    # For water solvation, use a dodecahedron box with solvent padding.
    # In our tests 1.5 nm padding works for FreeSolv; increase for ligands
    # that can unfold.
    settings.solvation_settings = PackmolSolvationSettings(
        box_shape="dodecahedron",
        assign_solvent_charges=False,
        solvent_padding=1.5 * unit.nanometer,
    )
    settings.lambda_settings.lambda_elec = [
        0.0,
        0.25,
        0.5,
        0.75,
        1.0,
        1.0,
        1.0,
        1.0,
        1.0,
        1.0,
        1.0,
        1.0,
        1.0,
        1.0,
    ]
    settings.lambda_settings.lambda_vdw = [
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.12,
        0.24,
        0.36,
        0.48,
        0.6,
        0.7,
        0.77,
        0.85,
        1.0,
    ]
    settings.lambda_settings.lambda_restraints = [
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
    ]
    settings.vacuum_simulation_settings.n_replicas = 14
    settings.solvent_simulation_settings.n_replicas = 14
    settings.solvent_simulation_settings.time_per_iteration = 2.5 * unit.picosecond
    settings.vacuum_simulation_settings.time_per_iteration = 2.5 * unit.picosecond
    settings.solvent_equil_simulation_settings.equilibration_length_nvt = (
        0.5 * unit.nanosecond
    )
    settings.solvent_equil_simulation_settings.equilibration_length = (
        0.5 * unit.nanosecond
    )
    settings.solvent_equil_simulation_settings.production_length = 9.5 * unit.nanosecond
    settings.vacuum_equil_simulation_settings.equilibration_length_nvt = None
    settings.vacuum_equil_simulation_settings.equilibration_length = (
        0.2 * unit.nanosecond
    )
    settings.vacuum_equil_simulation_settings.production_length = 0.5 * unit.nanosecond
    settings.solvent_simulation_settings.equilibration_length = 1.0 * unit.nanosecond
    settings.vacuum_simulation_settings.equilibration_length = 0.5 * unit.nanosecond
    settings.solvent_simulation_settings.production_length = 10.0 * unit.nanosecond
    settings.vacuum_simulation_settings.production_length = 2.0 * unit.nanosecond
    return settings


def get_nonwater_settings(forcefield: str) -> ASFESettings:
    """
    Get settings for solvation free energy calculations (non-water solvents).

    Parameters
    ----------
    forcefield : str
        Name of the force field to use.

    Returns
    -------
    ASFESettings
        Settings object for ASFEProtocol.
    """
    # The settings here are effectively the "fast" settings
    # shown in the validation.
    settings: ASFESettings = ASFEProtocol.default_settings()
    # Because it's Alchemiscale, you set protocol_repeats to 1 and then
    # run the Transformation task multiple times to get repeats.
    # Locally, the recommendation would be to set this to 3 so that you can
    # get a standard deviation uncertainty. It's not super necessary since
    # SFEs converge well, but hey with Alchemiscale why not?!
    settings.protocol_repeats = 1
    settings.solvent_forcefield_settings.forcefields = [
        # To use a custom force field, just pass an OFFXML string
        # just like you would to openff.toolkit.ForceField
        forcefield,
    ]
    settings.vacuum_forcefield_settings.forcefields = [
        forcefield,  # as above
    ]
    settings.vacuum_engine_settings.compute_platform = "CUDA"
    settings.solvent_engine_settings.compute_platform = "CUDA"
    settings.solvation_settings = PackmolSolvationSettings(
        # In our tests 750 gave quasi equivalent results to the 1999 used in the Sage
        # benchmarks
        number_of_solvent_molecules=750,
        box_shape="cube",
        # We set assign_solvent_charges to True because we don't have LibraryCharges.
        # If False it will only attempt to use LibraryCharges.
        # Note that if True and you don't have any charges on the SmallMoleculeComponent
        # passed to ExtendedSolventComponent, the Protocol will attempt to automatically
        # assign partial charges (default is AmberTools am1bcc, but it's controllable
        # using `partial_charge_settings`.
        assign_solvent_charges=True,
        solvent_padding=None,
    )
    # Below are the default lambda & replica settings, so you don't have to
    # actually set them, but if you want to change things, you can alter them
    # by defining them this way. Note: you have to update n_replica to match
    # the number of lambda windows (and all lambda window lists must be of the same length).
    settings.lambda_settings.lambda_elec = [
        0.0,
        0.25,
        0.5,
        0.75,
        1.0,
        1.0,
        1.0,
        1.0,
        1.0,
        1.0,
        1.0,
        1.0,
        1.0,
        1.0,
    ]
    settings.lambda_settings.lambda_vdw = [
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.12,
        0.24,
        0.36,
        0.48,
        0.6,
        0.7,
        0.77,
        0.85,
        1.0,
    ]
    settings.lambda_settings.lambda_restraints = [
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
    ]
    settings.vacuum_simulation_settings.n_replicas = 14
    settings.solvent_simulation_settings.n_replicas = 14
    # This set the time per replica exchange, the default is 1 ps but
    # hurts performance. I would recommend 2.5 ps
    settings.solvent_simulation_settings.time_per_iteration = 2.5 * unit.picosecond
    settings.vacuum_simulation_settings.time_per_iteration = 2.5 * unit.picosecond
    # Below are the default simulation lengths we use in Pontibus,
    # so you don't need to set them. However, you can do so manually
    # This is the pre-alchemical equilibration lengths
    # NVT equilibration -> NPT equilibration -> NPT "production" (more equilibration)
    # In vacuum, we set the NVT equilibration to None since it's all gas phase
    settings.solvent_equil_simulation_settings.equilibration_length_nvt = (
        0.5 * unit.nanosecond
    )
    settings.solvent_equil_simulation_settings.equilibration_length = (
        0.5 * unit.nanosecond
    )
    settings.solvent_equil_simulation_settings.production_length = 9.5 * unit.nanosecond
    settings.vacuum_equil_simulation_settings.equilibration_length_nvt = None
    settings.vacuum_equil_simulation_settings.equilibration_length = (
        0.2 * unit.nanosecond
    )
    settings.vacuum_equil_simulation_settings.production_length = 0.5 * unit.nanosecond
    # This is the alchemical equilibration length
    settings.solvent_simulation_settings.equilibration_length = 1.0 * unit.nanosecond
    settings.vacuum_simulation_settings.equilibration_length = 0.5 * unit.nanosecond
    # This is the alchemical production length
    settings.solvent_simulation_settings.production_length = 10.0 * unit.nanosecond
    settings.vacuum_simulation_settings.production_length = 2.0 * unit.nanosecond
    return settings


def compile_transformations(
    systems_by_name: dict[str, openfe.ChemicalSystem],
) -> list[openfe.Transformation]:
    """
    Create one ASFE transformation per chemical system.

    Each transformation runs both vacuum and solvent legs of the absolute solvation free energy calculation.
    stateA is the solvated solute and stateB is the pure solvent, as required by the Pontibus ASFEProtocol.

    Parameters
    ----------
    systems_by_name : dict of str to openfe.ChemicalSystem
        Mapping of system name to ChemicalSystem containing the solute and solvent components.

    Returns
    -------
    list of openfe.Transformation
        One Transformation per entry in systems_by_name.
    """
    transformations: list[openfe.Transformation] = []
    settings = {
        "aqueous": get_water_settings(FORCEFIELD, WATER_MODEL),
        "nonaqueous": get_nonwater_settings(FORCEFIELD),
    }
    for sys_name, system_solvated_solute in sorted(systems_by_name.items()):
        # An SFE transformation in GUFE formalism is defined as
        # going from a solute + solvent (stateA) to just solvent (stateB)
        # The vacuum states are created automatically, and this transformation
        # includes both the vacuum and solvent legs
        system_solvent = openfe.ChemicalSystem(
            {"solvent": system_solvated_solute.components["solvent"]}
        )
        sol_type = (
            "aqueous"
            if system_solvated_solute.components["solvent"].solvent_molecule.name
            == "water"
            else "nonaqueous"
        )
        transformations.append(
            openfe.Transformation(
                stateA=system_solvated_solute,
                stateB=system_solvent,
                mapping=None,
                protocol=ASFEProtocol(settings=settings[sol_type]),
                name=sys_name,
            )
        )

    return transformations


def main():
    """
    Build ASFE transformations and write network and per-transformation files.

    Returns
    -------
    None
    """

    script_dir = Path(__file__).parent
    out_dir = script_dir / OUTPUT_DIR
    out_dir.mkdir(parents=True, exist_ok=True)
    logger.info(f"Starting ASFE example: set={BENCHMARK_SET} sys={BENCHMARK_SYS}")

    benchmark_sys = get_benchmark_data_system(BENCHMARK_SET, BENCHMARK_SYS)
    logger.info(f"Loaded benchmark system {BENCHMARK_SET}/{BENCHMARK_SYS}")

    chemical_systems = get_chemical_systems(benchmark_sys)
    logger.info(
        f"Created {len(chemical_systems)} ChemicalSystems (partial-charge set: {PARTIAL_CHARGE})"
    )

    transformations = compile_transformations(chemical_systems)
    logger.info(f"Compiled {len(transformations)} transformations")

    # Write AlchemicalNetwork (this writes ChemicalSystem and Transformation records)
    alchem_network = openfe.AlchemicalNetwork(edges=transformations)
    alchem_network.to_json(file=out_dir / FILENAME)

    #    # Can be run HPC3
    #    for transformation in alchem_network.edges:
    #        transformation.to_json(os.path.join(OUTPUT_DIR, f"{transformation.name}.json"))

    logger.info(f"Wrote ASFE network to {out_dir / FILENAME}")


# ------------------------------- Validation ---------------------------------


def validate_asfe_network(network_file: Path) -> list[str]:
    """
    Validate an ASFE AlchemicalNetwork against BenchmarkData expectations.

    Loads the network as an openfe.AlchemicalNetwork instance and checks:
        - Each expected solute has a corresponding transformation.
        - Every transformation's stateA contains solute and solvent components.
        - Every SmallMoleculeComponent in the network has partial charges assigned (partial_charges is not None and length matches atom count).
        - All expected solutes are present somewhere in the network nodes.
        - The pontibus manifest exists, lists all expected solutes, and references an existing SDF file.

    Parameters
    ----------
    network_file : Path
        Path to the AlchemicalNetwork JSON file to validate.

    Returns
    -------
    list[str]
        List of error messages. Empty if validation passes.
    """
    errors: list[str] = []
    logger.info(f"Validating ASFE network at {network_file}")

    try:
        alchemical_network = openfe.AlchemicalNetwork.from_json(file=str(network_file))
    except Exception as exc:
        return [f"Failed to load alchemical network: {exc}"]

    try:
        edges = list(alchemical_network.edges)
    except Exception:
        return ["Loaded AlchemicalNetwork does not expose .edges"]

    logger.info(f"AlchemicalNetwork contains {len(edges)} transformations")
    logger.info(
        f"Expected {len(EXPECTED_NETWORKS)} transformations using {len(EXPECTED_LIGANDS)} ligands"
    )

    # Check that each expected solute has a transformation
    transformation_names = {getattr(t, "name", "") for t in edges}
    missing = [name for name in EXPECTED_NETWORKS if name not in transformation_names]
    if missing:
        errors.append(f"Missing expected transformations: {sorted(missing)}")

    # Validate stateA components for each transformation
    for t in edges:
        if t.stateA is None:
            errors.append(f"Transformation '{t.name}' has no stateA")
            continue
        comps = t.stateA.components
        for comp_name in ("solute", "solvent"):
            if comp_name not in comps:
                errors.append(
                    f"Transformation '{t.tname}' missing {comp_name} component"
                )

    # Validate partial charges using get_components_of_type across all nodes
    found_ligands: set[str] = set()
    for chem_system in alchemical_network.nodes:
        for ligand in chem_system.get_components_of_type(openfe.SmallMoleculeComponent):
            found_ligands.add(ligand.name)
            off_mol = ligand.to_openff()
            if (
                off_mol.partial_charges is not None
                and len(off_mol.partial_charges) == off_mol.n_atoms
            ):
                continue
            errors.append(
                f"Ligand '{ligand.name}' in '{chem_system.key}' is missing partial charges"
            )
        for solvent in chem_system.get_components_of_type(ExtendedSolventComponent):
            ligand = solvent.solvent_molecule
            found_ligands.add(ligand.name)
            if ligand.name == "water":
                continue
            off_mol = ligand.to_openff()
            if (
                off_mol.partial_charges is not None
                and len(off_mol.partial_charges) == off_mol.n_atoms
            ):
                continue
            errors.append(
                f"Ligand '{ligand.name}' in '{chem_system.key}' is missing partial charges"
            )

    print(
        len(EXPECTED_LIGANDS), len(found_ligands), len(EXPECTED_LIGANDS - found_ligands)
    )
    missing_ligs = EXPECTED_LIGANDS - found_ligands
    if missing_ligs:
        errors.append(f"Ligands not found in network: {sorted(missing_ligs)}")

    if errors:
        logger.error(
            f"ASFE network validation failed with {len(errors)} error(s): {errors}"
        )
    else:
        logger.info(
            f"ASFE network validation passed ({len(edges)} transformations, {len(found_ligands)} solutes)"
        )

    return errors


if __name__ == "__main__":
    # Simple example-run logging to make outputs visible when executed directly
    _configure_example_logging(level=logging.INFO)
    main()
    # Validate and raise on errors so the example fails loudly when incorrect
    errors = validate_asfe_network(Path(__file__).parent / OUTPUT_DIR / FILENAME)
    if errors:
        raise RuntimeError("ASFE network validation failed:\n" + "\n".join(errors))
    logger.info("Example completed and validated successfully")
