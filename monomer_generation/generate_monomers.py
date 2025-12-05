"""
Generates jsons files using the new format up-to-date as of 3/14/23
"""

import enum
import json
import multiprocessing
import os
import signal
import sys
import time
from collections.abc import Generator
from pathlib import Path
from types import FrameType
import logging

import openmm
import tqdm
import tqdm.contrib
from tqdm.contrib.logging import logging_redirect_tqdm
from monomer_smiles_input import ALL_SMILES_INPUT
from openff.toolkit import Topology
from openff.toolkit.utils.toolkits import (
    GLOBAL_TOOLKIT_REGISTRY,
    OpenEyeToolkitWrapper,
)
from openff.units.openmm import to_openmm as to_openmm_quantity
from simtk import unit
from substructure_generator import SubstructureGenerator

sys.path.append(os.path.abspath(__file__ + "/../.."))  # TODO: fix this mess
from pdb_file_search import PDBFiles


class Columns(enum.IntEnum):
    FILENAME = 0
    SUCCESS = 1
    EXCEPTION_MESSAGE = 2
    TIME = 3


# CONFIGURATION CONSTANTS
WRITE_RESULTING_TOPOLOGIES: bool = False
SKIP_EXISTING_JSONTOPS: bool = False
N_PROCESSES: int | None = 8  # None: Auto-detect
SORT_BY: tuple[Columns, ...] = (Columns.SUCCESS, Columns.TIME)
PDB_FILE_TIMEOUT_SEC: None | int = None  # None: No timeout.
# END OF CONFIGURATION CONSTANTS

# Note that PDB_FILE_TIMEOUT_SEC uses signal.alarm, which can't interrupt in the
# middle of a Python bytecode instruction. This means that if the process is
# in a long-running C calculation (such as matching SMARTS to an RDMol, or
# splitting up an RDMol into its constituent molecules), the timeout will
# fire at the end of that calculation. If the process is taking longer than you
# expect, you can use ctrl+c to interrupt - you'll still get results for what's
# already finished

GLOBAL_TOOLKIT_REGISTRY.deregister_toolkit(OpenEyeToolkitWrapper)


def timeout_handler(signum: int, frame: FrameType | None):
    raise TimeoutError(f"Timeout after {PDB_FILE_TIMEOUT_SEC} seconds")


def _identify_all_molecules(
    self,
) -> dict[int, tuple[int, dict[int, int]]]:
    identity_maps: dict[int, tuple[int, dict[int, int]]] = dict()
    already_matched_mols = set()

    for mol1_idx in range(self.n_molecules):
        if mol1_idx in already_matched_mols:
            continue
        mol1 = self.molecule(mol1_idx)
        identity_maps[mol1_idx] = (
            mol1_idx,
            {i: i for i in range(mol1.n_atoms)},
        )

    return identity_maps


def minimize_energy(off_topology, forcefield, max_iters):
    # mol should already have one conformer...

    # pdbfile = PDBFile(pdbfile)
    # omm_topology = pdbfile.topology
    omm_topology = off_topology.to_openmm()

    for m in off_topology.molecules:
        m.assign_partial_charges(partial_charge_method="gasteiger")

    start = time.time()
    system = forcefield.create_openmm_system(
        off_topology,
        allow_nonintegral_charges=True,
        charge_from_molecules=off_topology.molecules,
    )
    time_to_parameterize = time.time() - start

    time_step = 2 * unit.femtoseconds  # simulation timestep
    temperature = 1000 * unit.kelvin  # simulation temperature
    friction = 1 / unit.picosecond  # collision rate
    integrator = openmm.LangevinIntegrator(temperature, friction, time_step)
    simulation = openmm.app.Simulation(omm_topology, system, integrator)
    openmm_positions = to_openmm_quantity(off_topology.get_positions())
    simulation.context.setPositions(openmm_positions)

    start = time.time()
    simulation.minimizeEnergy(maxIterations=max_iters)
    time_to_energy_minimize = time.time() - start
    simulation.step(2)
    st = simulation.context.getState(getPositions=True, getEnergy=True)
    return st.getPotentialEnergy(), time_to_parameterize, time_to_energy_minimize


def load_file(
    file_name: str,
    monomer_info: dict[str, tuple[str, list[int]]],
    test_load: bool = True,
) -> Generator[tuple[Path, bool, str | None, float]]:
    # time_to_load
    # time_to_parameterize
    # time_to_energy_minimize
    engine = SubstructureGenerator()
    for name, substructure_and_caps in monomer_info.items():
        smarts, caps = substructure_and_caps
        if caps:
            engine.add_monomer_as_smarts_fragment(smarts, name, caps)
        else:
            engine.add_monomer(name, smarts)
    json_file = json_dir / Path(file_name + ".json")
    engine.output_monomer_info_json(json_file)

    pdb_files = PDBFiles.search(file_name)
    for pdb_file in pdb_files:
        if test_load:
            assert pdb_file != None
            substructs = engine.get_monomer_info_dict()["monomers"]

            # manually ensure that no molecules are cached to obtain the worst-case time complexity
            # for if parameterization should be done over ALL atoms, since that is the time-complexity
            # problem we are interested in solving. This can be done with some manipulation of
            # openff's identical_molecule_groups property for version 0.13.2
            Topology._identify_chemically_identical_molecules = _identify_all_molecules

            current_dir = Path(__file__).parent.parent.resolve()
            pdb_out = (
                current_dir
                / "results"
                / Path(pdb_file).relative_to(current_dir / "compatible_pdbs")
            )
            monomers_out = pdb_out.with_suffix(".monomers.json")
            jsontop_out = pdb_out.with_suffix(".topology.json")
            pdb_file_short = pdb_file.relative_to(pdb_file.parent.parent)

            if WRITE_RESULTING_TOPOLOGIES:
                pdb_out.parent.mkdir(parents=True, exist_ok=True)
                monomers_out.write_text(json.dumps(substructs))
                pdb_out.write_text(pdb_file.read_text())

            if SKIP_EXISTING_JSONTOPS and jsontop_out.exists():
                continue

            start = time.time()
            if PDB_FILE_TIMEOUT_SEC is not None:
                signal.signal(signal.SIGALRM, timeout_handler)
                signal.alarm(PDB_FILE_TIMEOUT_SEC)
            try:
                top: Topology = Topology.from_pdb(
                    str(pdb_file),
                    _custom_substructures=substructs,
                )
            except Exception as e:
                success = False
                exc = str(e)
            else:
                success = successfully_loaded(top)
                exc = None
                if WRITE_RESULTING_TOPOLOGIES and success:
                    jsontop_out.write_text(top.to_json())
            finally:
                if PDB_FILE_TIMEOUT_SEC is not None:
                    signal.alarm(0)
            time_to_load = time.time() - start

            yield (pdb_file_short, success, exc, time_to_load)

            # # desolvate since not all systems have solvent
            # new_top = Topology()
            # for mol in top.molecules:
            #     if mol != Molecule.from_smiles("[H]-O-[H]"):
            #         new_top.add_molecule(mol)
            # top = new_top

            # num_atoms = top.n_atoms

            # general_offxml = 'openff-2.0.0.offxml'
            # amber_offxml = 'ff14sb_off_impropers_0.0.3.offxml'
            # water_model = 'tip3p_fb-1.1.0.offxml'
            # forcefield = ForceField(general_offxml, amber_offxml, water_model)
            # # forcefield = ForceField(general_offxml, water_model)
            # forcefield.deregister_parameter_handler('ToolkitAM1BCC')
            # forcefield.get_parameter_handler('ChargeIncrementModel', {"version":0.3, "partial_charge_method":"gasteiger"})

            # energy = np.nan
            # for max_iters in [0]:
            #     try:
            #         print(f"trying energy minimization with {max_iters} maximum iterations")
            #         energy, time_to_parameterize, time_to_energy_minimize = minimize_energy(top, forcefield, max_iters)
            #         break # if successfully minimized without coordinate explosion
            #     except openmm.OpenMMException as e:
            #         print(f"openmm exception: {e}")

            # print(energy)
            # # with open("polymer_energies.txt", "a") as file:
            # #     file.write(f"{pdb_file.stem}, {num_atoms}, {energy}, {time_to_load}, {time_to_parameterize}, {time_to_energy_minimize}\n")
            # # #______________________________________________________________________________


def load_file_star(
    tup: tuple[str, dict[str, tuple[str, list[int]]]],
) -> list[tuple[Path, bool, str | None, float]]:
    return list(load_file(*tup))


def successfully_loaded(top: Topology) -> bool:
    match_info = [atom.metadata["match_info"] for atom in top.atoms]
    return all([bool(match) for match in match_info])


if __name__ == "__main__":
    # Make a file to store new jsons (TODO: change this to any new file structure)
    current_dir = Path(__file__).parent.resolve()
    json_dir = current_dir / Path("json_files")
    json_dir.mkdir(parents=False, exist_ok=True)
    os.chdir(current_dir)

    # create object for json creation and loading:
    # with open("polymer_energies.txt", "w") as file:
    #     file.write("name, num_atoms, energy, time_to_load, time_to_parameterize, time_to_energy_minimize\n")

    results: list[tuple[Path, bool, str | None, float]] = []
    interrupt = False
    with multiprocessing.Pool(N_PROCESSES) as pool:
        try:
            for elems in tqdm.tqdm(
                pool.imap_unordered(
                    load_file_star,
                    [*ALL_SMILES_INPUT.items()],
                ),
                total=len(ALL_SMILES_INPUT),
            ):
                results.extend(elems)
        except KeyboardInterrupt:
            interrupt = True

    longest_pdb_file = max(len(str(pdb_file)) for pdb_file, *_ in results)

    n_successes = 0
    n_failures = 0
    for pdb_file, success, exc, time_to_load in sorted(
        results,
        key=lambda tup: tuple(tup[i] for i in SORT_BY),
    ):
        if success:
            msg = ""
        elif exc is None:
            msg = ": topology loaded but failed to validate"
        else:
            msg = f": {exc}"

        if success:
            n_successes += 1
            success_str = "SUCCESS"
        else:
            n_failures += 1
            success_str = "FAILURE"

        print(
            f"{success_str}: {str(pdb_file):{longest_pdb_file}} in {time_to_load:10.3f} s{msg}",
        )

    print(f"There were {n_successes} successes out of {len(results)} completed PDB files!")
    if interrupt:
        print("(processing was interrupted by user)")
