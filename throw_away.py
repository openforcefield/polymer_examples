from openff.toolkit import Topology
from openff.toolkit.utils import get_data_file_path
from rdkit import Chem

middle_monomer = "*CC*" # simple smiles repr of monomer 
end_group = "*CC" # "*" represents neighboring atoms

templates = []
for monomer in [middle_monomer, end_group]:
    rdmol = Chem.MolFromSmiles(monomer) # load as a molecule to add Hs
    rdmol = Chem.AddHs(rdmol)           # Add explicit Hs
    rdmol = Chem.MolFromSmarts(Chem.MolToSmarts(rdmol)) # reload as query
    for atom in rdmol.GetAtoms():
        if atom.GetAtomicNum() > 0: # internal atoms must have #<>, D<>, and +<>
            a_num = atom.GetAtomicNum()
            D_num = len([0 for _ in atom.GetBonds()]) # explicit bonds 
            F_num = atom.GetFormalCharge()
            query_string = f"[#{a_num}D{D_num}{F_num:+}:{atom.GetAtomMapNum()}]"
            query = Chem.AtomFromSmarts(query_string) # manually set query
            atom.SetQuery(query)
        else: # Atomic Number = 0, which indicates a neighbor atom
            query_string = f"[*:{atom.GetAtomMapNum()}]"
            query = Chem.AtomFromSmarts(query_string)
            atom.SetQuery(query)

    # Finally, set atom map numbers
    [atom.SetAtomMapNum(atom.GetIdx() + 1) for atom in rdmol.GetAtoms()]

    # output the template and remove redundant "&" chars
    smarts_string = Chem.MolToSmarts(rdmol)
    smarts_string = smarts_string.replace('&', '')
    templates.append(smarts_string)

print("Templates generated as:")
[print(temp) for temp in templates]

pdb_file =  get_data_file_path("systems/test_systems/PE.pdb")
top = Topology.from_pdb(pdb_file, _custom_substructures={"PE": templates})
print("Topology Successfully Loaded")