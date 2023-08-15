"""
Generates jsons files using the new format up-to-date as of 3/14/23
"""

from openff.toolkit import Topology, Molecule
import sys
import os
import numpy as np
from pathlib import Path
from rdkit import Chem
import time
import openmm
from copy import deepcopy
from simtk import unit
from openff.toolkit.typing.engines.smirnoff import ForceField
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

def parameterize(off_topology, forcefield):
    # mol should already have one conformer...

    # pdbfile = PDBFile(pdbfile)
    # omm_topology = pdbfile.topology
    omm_topology = off_topology.to_openmm()

    start = time.time()
    system = forcefield.create_openmm_system(off_topology, allow_nonintegral_charges=True)
    time_to_parameterize = time.time() - start

    return time_to_parameterize

def test_load(name, top, out_file):
    
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
        time_to_parameterize = parameterize(top, forcefield)
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

with open("tenk_timed_tests.txt", "w") as file:
    file.write("name, num_atoms, time_to_parameterize\n")

# for onek_file in Path("onek_polymers").iterdir():
#     tenk_file = ""
#     for file in Path("tenk_polymers").iterdir():
#         if file.stem.split("_10000")[0] == onek_file.stem:
#             tenk_file = file
#     if not tenk_file:
#         print(f"failed to find both onek and tenk files for {onek_file.stem}")
#         continue

    # monomer_info = {}
    # for file_name, monomer_info_query in ALL_SMILES_INPUT.items():
    #     if file_name.split("_modified")[0] == onek_file.stem.split("_modified")[0]:
    #         monomer_info = monomer_info_query
    #         break
    # if not monomer_info:
    #     print(f"failed to find monomer info for {onek_file}")
    #     continue

    # engine = SubstructureGenerator()
    # for name, substructure_and_caps in monomer_info.items():
    #     smarts, caps = substructure_and_caps
    #     if caps:
    #         engine.add_monomer_as_smarts_fragment(smarts, name, caps)
    #     else:
    #         engine.add_monomer(name, smarts)

    # substructs = engine.get_monomer_info_dict()["monomers"]

    # test_load(onek_file, substructs, "onek_timed_tests.txt")

systems = [
           (2500, 1),
           (1250, 2), 
           (833, 3), 
           (625, 4),
           (500, 5)
          ]
for num_atoms, num_chains in systems:
    mol = Molecule.from_smiles("CC"*num_atoms)
    for i in range(0,3):
        # manually ensure that no molecules are cached to obtain the worst-case time complexity 
        # for if parameterization should be done over ALL atoms, since that is the time-complexity
        # problem we are interested in solving. This can be done with some manipulation of
        # openff's identical_molecule_groups property for version 0.13.2
        Topology._identify_chemically_identical_molecules = _identify_all_molecules
        top = Topology.from_molecules([mol]*num_chains)
        test_load(f"PE_{num_atoms}x{num_chains}", top, "tenk_timed_tests.txt")