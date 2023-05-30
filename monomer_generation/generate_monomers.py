"""
Generates jsons files using the new format up-to-date as of 3/14/23
"""

from openff.toolkit import Topology
from substructure_generator import SubstructureGenerator
import sys
import os
from monomer_smiles_input import ALL_SMILES_INPUT
from pathlib import Path

sys.path.append(os.path.abspath(__file__ + "/../..")) # TODO: fix this mess
from pdb_file_search import PDBFiles

def successfully_loaded(top):
    match_info = [atom.metadata["match_info"] for atom in top.atoms]
    return all([bool(match) for match in match_info])

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
    if test_load:
        assert pdb_file != None
        substructs = engine.get_monomer_info_dict()["monomers"]
        top = Topology.from_pdb(str(pdb_file), _custom_substructures=substructs)
        assert successfully_loaded(top)