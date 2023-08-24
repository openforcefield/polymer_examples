"""
Generates jsons files using the new format up-to-date as of 3/14/23
"""

from openff.toolkit import Topology, Molecule
from substructure_generator import SubstructureGenerator
import sys
import os
import numpy as np
from monomer_smiles_input import ALL_SMILES_INPUT
from pathlib import Path
from partition import partition
from rdkit import Chem
from openeye import oechem
import time
import openmm
from openmm.app import PDBFile
from copy import deepcopy
from simtk import unit
from openff.toolkit.typing.engines.smirnoff import ForceField
from openff.toolkit.typing.engines.smirnoff import vdWHandler, AngleHandler
from openff.toolkit.utils.toolkits import GLOBAL_TOOLKIT_REGISTRY, ToolkitRegistry, OpenEyeToolkitWrapper
from openff.units.openmm import to_openmm as to_openmm_quantity

from typing import (
    Dict,
    List,
    Tuple,
)

GLOBAL_TOOLKIT_REGISTRY.deregister_toolkit(OpenEyeToolkitWrapper)

def _identify_all_molecules(
    self,
    ) -> Dict[int, Tuple[int, Dict[int, int]]]:
    identity_maps: Dict[int, Tuple[int, Dict[int, int]]] = dict()
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
    system = forcefield.create_openmm_system(off_topology, allow_nonintegral_charges=True, charge_from_molecules=off_topology.molecules)
    time_to_parameterize = time.time() - start

    time_step = 2*unit.femtoseconds  # simulation timestep
    temperature = 1000*unit.kelvin  # simulation temperature
    friction = 1/unit.picosecond  # collision rate
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

sys.path.append(os.path.abspath(__file__ + "/../..")) # TODO: fix this mess
from pdb_file_search import PDBFiles

def successfully_loaded(top):
    match_info = [atom.metadata["match_info"] for atom in top.atoms]
    return all([bool(match) for match in match_info])

# Make a file to store new jsons (TODO: change this to any new file structure)
current_dir = Path(__file__).parent.resolve()
json_dir = current_dir / Path("json_files")
json_dir.mkdir(parents=False, exist_ok=True)
os.chdir(current_dir)

# set flag if the script should try to test_load the new json file
test_load = True

# create object for json creation and loading:
# with open("polymer_energies.txt", "w") as file:
#     file.write("name, num_atoms, energy, time_to_load, time_to_parameterize, time_to_energy_minimize\n")

for file_name, monomer_info in ALL_SMILES_INPUT.items():
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
            print(f"testing {file_name}")
            assert pdb_file != None
            substructs = engine.get_monomer_info_dict()["monomers"]

            # manually ensure that no molecules are cached to obtain the worst-case time complexity 
            # for if parameterization should be done over ALL atoms, since that is the time-complexity
            # problem we are interested in solving. This can be done with some manipulation of
            # openff's identical_molecule_groups property for version 0.13.2
            Topology._identify_chemically_identical_molecules = _identify_all_molecules

            start = time.time()
            top = Topology.from_pdb(str(pdb_file), _custom_substructures=substructs)
            time_to_load = time.time() - start
            assert successfully_loaded(top)

            # desolvate since not all systems have solvent
            new_top = Topology()
            for mol in top.molecules:
                if mol != Molecule.from_smiles("[H]-O-[H]"):
                    new_top.add_molecule(mol)
            top = new_top

            num_atoms = top.n_atoms

            general_offxml = 'openff-2.0.0.offxml'
            amber_offxml = 'ff14sb_off_impropers_0.0.3.offxml'
            water_model = 'tip3p_fb-1.1.0.offxml'
            forcefield = ForceField(general_offxml, amber_offxml, water_model)
            # forcefield = ForceField(general_offxml, water_model)
            forcefield.deregister_parameter_handler('ToolkitAM1BCC')
            forcefield.get_parameter_handler('ChargeIncrementModel', {"version":0.3, "partial_charge_method":"gasteiger"})

            energy = np.nan
            for max_iters in [0]:
                try:
                    print(f"trying energy minimization with {max_iters} maximum iterations")
                    energy, time_to_parameterize, time_to_energy_minimize = minimize_energy(top, forcefield, max_iters)
                    break # if successfully minimized without coordinate explosion
                except openmm.OpenMMException as e:
                    print(f"openmm exception: {e}")

            print(energy)
            # with open("polymer_energies.txt", "a") as file:
            #     file.write(f"{pdb_file.stem}, {num_atoms}, {energy}, {time_to_load}, {time_to_parameterize}, {time_to_energy_minimize}\n")
            # #______________________________________________________________________________