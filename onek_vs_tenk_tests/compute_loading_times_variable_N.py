"""
Generates jsons files using the new format up-to-date as of 3/14/23
"""

from openff.toolkit import Topology, Molecule
import sys
import os
import numpy as np
from pathlib import Path
import pstats
from rdkit import Chem
import time
import openmm
from copy import deepcopy
from simtk import unit
from openff.toolkit.typing.engines.smirnoff import ForceField
from openff.toolkit.utils.toolkits import GLOBAL_TOOLKIT_REGISTRY, ToolkitRegistry, OpenEyeToolkitWrapper
from openff.units.openmm import to_openmm as to_openmm_quantity
import cProfile as profile

from typing import (
    Dict,
    List,
    Tuple,
)

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

def parameterize(off_topology, forcefield, prof):
    # mol should already have one conformer...

    # pdbfile = PDBFile(pdbfile)
    # omm_topology = pdbfile.topology
    omm_topology = off_topology.to_openmm()

    for m in off_topology.molecules:
        m.assign_partial_charges(partial_charge_method="gasteiger")

    start = time.time()
    prof.enable()
    system = forcefield.create_openmm_system(off_topology, allow_nonintegral_charges=True, charge_from_molecules=off_topology.molecules)
    prof.disable()
    time_to_parameterize = time.time() - start

    return time_to_parameterize

def test_load(name, top, out_file, prof):
    
    # try:
    #     start = time.time()
    #     top = Topology.from_pdb(str(pdb_file), _custom_substructures=substructs)
    #     time_to_load = time.time() - start
    # except Exception as e:
    #     print("exception during loading:")
    #     print(e)
    #     return 
    # assert successfully_loaded(top)

    # # desolvate since not all systems have solvent
    # new_top = Topology()
    # for mol in top.molecules:
    #     if mol != Molecule.from_smiles("[H]-O-[H]"):
    #         new_top.add_molecule(mol)
    # top = new_top

    num_atoms = top.n_atoms

    general_offxml = 'openff-2.0.0.offxml'
    amber_offxml = 'ff14sb_off_impropers_0.0.3.offxml'
    water_model = 'tip3p_fb-1.1.0.offxml'
    forcefield = ForceField(general_offxml, amber_offxml, water_model)
    # forcefield = ForceField(general_offxml, water_model)
    forcefield.deregister_parameter_handler('ToolkitAM1BCC')
    forcefield.get_parameter_handler('ChargeIncrementModel', {"version":0.3, "partial_charge_method":"gasteiger"})

    energy = np.nan
    try:
        time_to_parameterize = parameterize(top, forcefield, prof)
    except Exception as e:
        print(e)
        return

    print(energy)
    with open(out_file, "a") as file:
        file.write(f"{name}, {num_atoms}, {time_to_parameterize}\n")
    #______________________________________________________________________________

sys.path.append(os.path.abspath(__file__ + "/../.."))
sys.path.append(os.path.abspath(__file__ + "/../../monomer_generation"))

from pdb_file_search import PDBFiles
from substructure_generator import SubstructureGenerator
from monomer_smiles_input import ALL_SMILES_INPUT

def successfully_loaded(top):
    match_info = [atom.metadata["match_info"] for atom in top.atoms]
    return all([bool(match) for match in match_info])

# Make a file to store new jsons (TODO: change this to any new file structure)
current_dir = Path(__file__).parent.resolve()
os.chdir(current_dir)

with open("tenk_timed_tests_variable_N.txt", "w") as file:
    file.write("name, num_atoms, time_to_parameterize\n")

Ns = [3000]
for dop in Ns:
    for num_chains in [10]:
        dop_per_chain = int(dop / num_chains)
        mol = Molecule.from_smiles("CC"*dop_per_chain)
        for i in [3]: # range(0,3):
            # manually ensure that no molecules are cached to obtain the worst-case time complexity 
            # for if parameterization should be done over ALL atoms, since that is the time-complexity
            # problem we are interested in solving. This can be done with some manipulation of
            # openff's identical_molecule_groups property for version 0.13.2
            Topology._identify_chemically_identical_molecules = _identify_all_molecules
            top = Topology.from_molecules([mol]*num_chains)
            entry_name = f"PE_{dop_per_chain}x{num_chains}"
            run_string = f"test_load(\"{entry_name}\", top, \"tenk_timed_tests_variable_N.txt\")"
            
            prof = profile.Profile()
            test_load(entry_name, top, "tenk_timed_tests_variable_N.txt", prof)
            stream = open(f'profile_data/{entry_name}_t{i}.txt', 'w')
            stats = pstats.Stats(prof, stream=stream)
            stats.sort_stats('cumtime')
            stats.strip_dirs()
            stats.print_stats()