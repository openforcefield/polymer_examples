from pathlib import Path
import os
from monomer_smiles_input import ALL_SMILES_INPUT
from rdkit import Chem

current_dir = Path(__file__).parent.resolve()
os.chdir(current_dir)

bonds = {Chem.BondType.SINGLE: Chem.MolFromSmarts("*-*").GetBondWithIdx(0),
         Chem.BondType.DOUBLE: Chem.MolFromSmarts("*=*").GetBondWithIdx(0),
         Chem.BondType.TRIPLE: Chem.MolFromSmarts("*#*").GetBondWithIdx(0)}

for file_name, monomer_info in ALL_SMILES_INPUT.items():
    print(file_name)
    for name, substructure_and_caps in monomer_info.items():
        smarts, caps = substructure_and_caps
        rdmol = Chem.MolFromSmiles(smarts)
        rdmol = Chem.AddHs(rdmol)
        if Chem.BondType.AROMATIC in [b.GetBondType() for b in rdmol.GetBonds()]:
            Chem.Kekulize(rdmol)
            # for bond in rdmol.GetBonds():
            #     bond.SetQuery(bonds[bond.GetBondType()])
            print(f"{name}\t\t {Chem.MolToSmiles(rdmol, kekuleSmiles=True, allBondsExplicit=True, allHsExplicit=True)}")