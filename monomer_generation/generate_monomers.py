"""
Generates jsons files using the new format up-to-date as of 3/14/23
"""

from substructure_generator import SubstructureGenerator
from substructure_vizualizer import ChemistryEngine
import sys
import os
from monomer_smiles_input import ALL_SMILES_INPUT
from pathlib import Path

sys.path.append(os.path.abspath(__file__ + "/../..")) # TODO: fix this mess
from pdb_file_search import PDBFiles

# Make a file to store new jsons (TODO: change this to any new file structure)
current_dir = Path(__file__).parent.resolve()
json_dir = current_dir / Path("json_files")
json_dir.mkdir(parents=False, exist_ok=True)

# set flag if the script should try to test_load the new json file
test_load = True

# create object for json creation and loading:
for file_name, monomer_info in ALL_SMILES_INPUT.items():
    engine = SubstructureGenerator()
    for name, substructure_and_caps in monomer_info.items():
        smarts, caps = substructure_and_caps
        if caps:
            engine.add_monomer_as_smarts_fragment(smarts, name, caps)
        else:
            engine.add_monomer(name, smarts)
    json_file = json_dir / Path(file_name + ".json")
    engine.output_monomer_info_json(json_file)

    pdb_file = PDBFiles.search(file_name)
    if test_load and pdb_file != None:
        chem_engine = ChemistryEngine(str(pdb_file))
        assigned_atoms, assigned_bonds, unassigned_atoms, unassigned_bonds, _, _ = chem_engine.test_polymer_load(str(json_file), verbose=False)
        if len(unassigned_atoms) == 0 and len(unassigned_bonds) == 0:
            print(f"{file_name} successfully created with full coverage of the molecule")
        else:
            print(f"--> ***WARNING*** {file_name} missed {unassigned_atoms} atoms and {unassigned_bonds} during loading ***WARNING*** <--")


