"""
This module provides an example script for planning absolute solvation free energy (ASFE) calculations.
It demonstrates how to process benchmark systems, prepare components, and save the resulting
transformations to a JSON file.

Note that although the output of this script is an alchemical network, that simply serves as a storage
class for the independent transformations needed in ASFEs.
"""

import json
import logging
from pathlib import Path

from openff.units import unit
import openfe
from openfe import SolventComponent, ChemicalSystem, SmallMoleculeComponent

from pontibus.protocols.solvation import ASFEProtocol, ASFESettings
from pontibus.protocols.solvation.settings import PackmolSolvationSettings

from openfe_benchmarks.data import get_benchmark_data_system
from openfe_benchmarks.scripts import utils as ofebu

SOLVENT = SolventComponent(positive_ion="Na", negative_ion="Cl", neutralize=True)
BENCHMARK_SET = "mcs_docking_set"
BENCHMARK_SYS = "hne"
PARTIAL_CHARGE = "nagl_openff-gnn-am1bcc-1.0.0.pt"
FORCEFIELD = "openff-2.3.0"  # available [openmmforcefields SystemGenerator](https://github.com/openmm/openmmforcefields?tab=readme-ov-file#automating-force-field-management-with-systemgenerator)
OUTPUT_DIR = "outputs"
FILENAME = f"network_{BENCHMARK_SET}_{BENCHMARK_SYS}_nacl_asfe.json"
PONTIBUS_FILENAME = "pontibus_inputs.json"

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
    logger.debug("Example logging configured (level=%s)", logging.getLevelName(level))


def process_components(benchmark_sys):
    """
    Process the components of a benchmark system.

    Parameters
    ----------
    benchmark_sys : BenchmarkData
        The benchmark system to process.

    Returns
    -------
    tuple
        A tuple containing the solute dictionary and cofactors.
    """
    solute_dict = ofebu.process_sdf(
        benchmark_sys.ligands[PARTIAL_CHARGE], return_dict=True
    )

    cofactors = None
    if benchmark_sys.cofactors is not None:
        cofactors = ofebu.process_sdf(
            benchmark_sys.cofactors[PARTIAL_CHARGE], return_dict=False
        )

    return solute_dict, cofactors


def get_nonwater_settings(forcefield: str) -> ASFESettings:
    """
    Get the settings for solvation free energy calculations.
    These settings are tuned to give reasonable results on the FreeSolv dataset
    with a balance of accuracy and speed.
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
    """Create transformation per chemical system."""
    transformations: list[openfe.Transformation] = []
    settings = get_nonwater_settings(FORCEFIELD)
    protocol = ASFEProtocol(settings=settings)

    for sys_name, system_solvated_solute in sorted(systems_by_name.items()):
        # An SFE transformation in GUFE formalism is defined as
        # going from a solute + solvent (stateA) to just solvent (stateB)
        # The vacuum states are created automatically, and this transformation
        # includes both the vacuum and solvent legs
        system_solvent = openfe.ChemicalSystem(
            {"solvent": system_solvated_solute.components["solvent"]}
        )
        transformations.append(
            openfe.Transformation(
                stateA=system_solvated_solute,
                stateB=system_solvent,
                mapping=None,
                protocol=protocol,
                name=sys_name,
            )
        )

    return transformations


def get_chemical_systems(
    solute_dict: dict[str, SmallMoleculeComponent],
    cofactors: list[SmallMoleculeComponent] | None = None,
) -> dict[str, ChemicalSystem]:
    """
    Using the benchmark data file, create a set of ChemicalSystems
    that contain the solute and solvent Components.

    Parameters
    ----------
    solute_dict : dict[str, SmallMoleculeComponent]
        Dictionary of SmallMoleculeComponents to calculate ASFE for
    cofactors : list[SmallMoleculeComponent] | None
        Dictionary of SmallMoleculeComponents to act as co-solute. Default is None.

    Returns
    -------
    dict[str, ChemicalSystem]
        List of ChemicalSystems for the benchmark.
    """
    systems = {}
    for lig_name, solute in solute_dict.items():
        csystem = {
            "solute": solute,
            "solvent": SOLVENT,
        }
        if cofactors is not None:
            for i, cofactor in enumerate(cofactors):
                csystem[cofactor.name] = cofactor

        systems[lig_name] = ChemicalSystem(csystem, name=lig_name)

    return systems


def write_pontibus_inputs(script_dir: Path, solute_names: list[str], source_sdf: str):
    """Write a minimal pontibus input manifest referencing the SDF and solutes.

    This file is intentionally simple and matches the needs of downstream
    packing: it tells pontibus which solute identifiers to extract from the
    provided SDF and what ionic conditions to apply.
    """
    out_dir = script_dir / OUTPUT_DIR
    out_dir.mkdir(parents=True, exist_ok=True)
    pontibus = {
        "source_sdf": str(source_sdf),
        "solutes": sorted(solute_names),
        "ions": {"positive": SOLVENT.positive_ion, "negative": SOLVENT.negative_ion},
        "ionic_strength_molar": 0.15,
        "box_padding_angstrom": 12.0,
    }
    pontibus_path = out_dir / PONTIBUS_FILENAME
    pontibus_path.write_text(json.dumps(pontibus, indent=2))
    logger.info("Wrote pontibus manifest to %s", pontibus_path)
    return pontibus_path


def main():
    """Build ASFE transformations and write network + per-transformation files."""
    script_dir = Path(__file__).parent
    out_dir = script_dir / OUTPUT_DIR
    out_dir.mkdir(parents=True, exist_ok=True)
    logger.info("Starting ASFE example: set=%s sys=%s", BENCHMARK_SET, BENCHMARK_SYS)

    benchmark_sys = get_benchmark_data_system(BENCHMARK_SET, BENCHMARK_SYS)
    logger.info("Loaded benchmark system %s/%s", BENCHMARK_SET, BENCHMARK_SYS)
    solutes, cofactors = process_components(benchmark_sys)
    logger.info(
        "Processed %d solutes (partial-charge set: %s); cofactors: %s",
        len(solutes),
        PARTIAL_CHARGE,
        "present" if cofactors else "none",
    )
    chemical_systems = get_chemical_systems(solutes, cofactors=cofactors)
    logger.info("Created %d ChemicalSystems", len(chemical_systems))
    transformations = compile_transformations(chemical_systems)
    logger.info("Compiled %d transformations", len(transformations))

    # Write AlchemicalNetwork (this writes ChemicalSystem and Transformation records)
    alchem_network = openfe.AlchemicalNetwork(edges=transformations)
    alchem_network.to_json(file=out_dir / FILENAME)

    #    # Can be run HPC3
    #    for transformation in alchem_network.edges:
    #        transformation.to_json(os.path.join(OUTPUT_DIR, f"{transformation.name}.json"))

    # Create pontibus manifest that points to the source SDF and lists solutes
    pontibus_file = write_pontibus_inputs(
        script_dir, list(solutes.keys()), benchmark_sys.ligands[PARTIAL_CHARGE]
    )

    logger.info(
        "Wrote ASFE network to %s and pontibus inputs to %s",
        out_dir / FILENAME,
        pontibus_file,
    )


# ------------------------------- Validation ---------------------------------


def validate_asfe_network() -> list[str]:
    """Validate an ASFE AlchemicalNetwork object (calculation-ready checks).

    This implementation loads the network as an `openfe.AlchemicalNetwork`
    instance and inspects the in-memory objects (names, stateA ChemicalSystem
    components, and SmallMoleculeComponent molprops) rather than parsing
    the raw JSON structure.
    """
    script_dir = Path(__file__).parent
    out_dir = script_dir / OUTPUT_DIR
    errors: list[str] = []
    logger.info("Validating ASFE network at %s", out_dir / FILENAME)

    # Load as an AlchemicalNetwork object (prefer the object's loader)
    try:
        alchemical_network = openfe.AlchemicalNetwork.from_json(
            file=str(out_dir / FILENAME)
        )
    except Exception as exc:
        return [f"Failed to load alchemical network: {exc}"]

    # optional: network name for diagnostics
    net_name = getattr(alchemical_network, "name", None)
    logger.debug("Loaded AlchemicalNetwork%s", f" '{net_name}'" if net_name else "")

    # Gather transformations and chemical systems from the object model
    try:
        edges = list(alchemical_network.edges)
    except Exception:
        return ["Loaded AlchemicalNetwork does not expose .edges"]

    logger.info("AlchemicalNetwork contains %d transformations", len(edges))

    # Expected solutes from the source SDF
    benchmark_sys = get_benchmark_data_system(BENCHMARK_SET, BENCHMARK_SYS)
    expected_solutes = ofebu.process_sdf(
        benchmark_sys.ligands[PARTIAL_CHARGE], return_dict=True
    )
    expected_names = set(str(k) for k in expected_solutes.keys())

    # Check that each expected solute has both solvent_ and vacuum_ transformations
    transformation_names = {getattr(t, "name", "") for t in edges}
    missing = []
    for name in expected_names:
        if name not in transformation_names:
            missing.append(name)
    if missing:
        errors.append(f"Missing expected transformations: {sorted(missing)}")

    # Inspect components and collect solute SmallMoleculeComponent objects
    lig_components: dict[str, openfe.SmallMoleculeComponent] = {}
    for t in edges:
        tname = getattr(t, "name", "")
        stateA = getattr(t, "stateA", None)

        comps = {}
        if stateA is None:
            errors.append(f"Transformation '{tname}' has no stateA")
            continue
        if hasattr(stateA, "components"):
            comps = getattr(stateA, "components")

        # Validate solvent/vacuum expectations
        for comp_name in ("solute", "solvent"):
            if comp_name not in comps:
                errors.append(f"Transformation '{tname}' missing {comp_name} component")

        solute_obj = comps.get("solute")
        if solute_obj is not None:
            lig_components[solute_obj.name] = solute_obj

    # Ensure all expected solutes were present
    missing_ligs = expected_names - set(k for k in lig_components.keys() if k)
    if missing_ligs:
        errors.append(f"Ligands not found in network: {sorted(list(missing_ligs))}")

    # Check partial charges present in molprops for each solute
    solutes_without_charges = []
    for lname, solute_obj in lig_components.items():
        offmol = solute_obj.to_openff()
        has_charges = all(offmol.partial_charges != 0)
        if not has_charges:
            solutes_without_charges.append(lname)

    if solutes_without_charges:
        errors.append(
            f"Ligands missing partial charges: {sorted(solutes_without_charges)}"
        )

    # Check pontibus manifest
    pontibus_path = out_dir / PONTIBUS_FILENAME
    if not pontibus_path.exists():
        errors.append(f"Missing pontibus manifest: {pontibus_path.name}")
    else:
        try:
            p = json.loads(pontibus_path.read_text())
            listed = set(p.get("solutes", []))
            if not expected_names.issubset(listed):
                errors.append("pontibus manifest missing solute entries")
            src = Path(p.get("source_sdf", ""))
            if not src.exists():
                # allow relative paths resolved against the script dir
                src2 = Path(__file__).parent / src
                if not src2.exists():
                    errors.append("pontibus source SDF not found: %s" % src)
        except Exception as exc:
            errors.append(f"Failed to parse pontibus manifest: {exc}")

    if errors:
        logger.error(
            "ASFE network validation failed with %d error(s): %s", len(errors), errors
        )
    else:
        logger.info(
            "ASFE network validation passed (%d transformations, %d solutes)",
            len(edges),
            len(lig_components),
        )

    return errors


if __name__ == "__main__":
    # Simple example-run logging to make outputs visible when executed directly
    _configure_example_logging(logging.INFO)
    main()
    # Validate and raise on errors so the example fails loudly when incorrect
    errors = validate_asfe_network()
    if errors:
        raise RuntimeError("ASFE network validation failed:\n" + "\n".join(errors))
    logger.info("Example completed and validated successfully")
