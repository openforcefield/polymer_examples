
from openff.toolkit import Topology, Molecule
from substructure_generator import SubstructureGenerator
import sys
import json
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
from collections import defaultdict

from typing import (
    Dict,
    List,
    Tuple,
)

def substructures_disagree(s1, s2):
    
    def _fuzzy_query(query):
        """return a copy of Query which is less specific:
        - ignore aromaticity and hybridization of atoms (i.e. [#6] not C)
        - ignore bond orders
        - ignore formal charges
        """
        from rdkit import Chem

        # it's tricky from the Python API to properly edit queries,
        # but you can do SetQuery on Atoms/Bonds to edit them quite powerfully
        generic = Chem.MolFromSmarts("**")
        generic_bond = generic.GetBondWithIdx(0)

        fuzzy = Chem.Mol(query)
        neighbor_idxs = []
        for idx, a in enumerate(fuzzy.GetAtoms()):
            a.SetFormalCharge(0)
            if a.GetAtomicNum() > 0:
                a.SetQuery(
                    Chem.AtomFromSmarts(f"[#{a.GetAtomicNum()}D{a.GetDegree()}]")
                )
            else:
                a.SetQuery(generic.GetAtomWithIdx(0))
            a.SetNoImplicit(True)
            if a.GetAtomicNum() == 0:
                neighbor_idxs.append(idx)
        for b in fuzzy.GetBonds():
            b.SetIsAromatic(False)
            b.SetBondType(Chem.rdchem.BondType.SINGLE)
            b.SetQuery(generic_bond)
        return fuzzy, neighbor_idxs
    
    def _get_symmetrical_groups(fuzzy_query, substruct):
        """Returns those atoms and bonds whose chemical information
        is ambiguous due to resonance forms or symmetrical groups. Conflicts
        in assignment are ignored for these atoms when two queries have the same
        atoms in resonance/symmetry"""
        from copy import deepcopy

        from rdkit import Chem

        qmol = deepcopy(fuzzy_query)
        for atom in qmol.GetAtoms():  # reset queries and map numbers
            atom.SetAtomMapNum(
                atom.GetIdx()
            )  # reorder atom map nums to later recover ids

        qmol = Chem.RemoveAllHs(qmol)
        idx_to_map_num = dict(
            [(a.GetIdx(), a.GetAtomMapNum()) for a in qmol.GetAtoms()]
        )
        automorphs = fuzzy_query.GetSubstructMatches(qmol, uniquify=0)
        ambiguous_bonds = []
        ambiguous_atoms = []
        for automorph in automorphs:
            # check for conflicting chemical information
            automorph = dict(
                [
                    (idx_to_map_num[idx], a)
                    for idx, a in enumerate(list(automorph))
                    if idx_to_map_num[idx] != a
                ]
            )  # only care about cases of different matching

            for atom_iso, new_atom_iso in automorph.items():
                atom = substruct.GetAtomWithIdx(atom_iso)
                new_atom = substruct.GetAtomWithIdx(new_atom_iso)
                # new_atom = substruct.GetAtomWithIdx(automorph[atom.GetIdx()])
                if atom.GetFormalCharge() != new_atom.GetFormalCharge():
                    if atom.GetIdx() not in ambiguous_atoms:
                        ambiguous_atoms.append(atom.GetIdx())

            for bond in substruct.GetBonds():
                if (
                    bond.GetBeginAtom().GetAtomicNum() == 1
                    or bond.GetEndAtom().GetAtomicNum() == 1
                ):  # we remove Hs for matching so must remove here as well
                    continue
                if (
                    bond.GetBeginAtomIdx() in automorph
                    or bond.GetEndAtomIdx() in automorph
                ):
                    new_bond_begin_idx = automorph.get(
                        bond.GetBeginAtomIdx(), bond.GetBeginAtomIdx()
                    )
                    new_bond_end_idx = automorph.get(
                        bond.GetEndAtomIdx(), bond.GetEndAtomIdx()
                    )
                    new_bond = substruct.GetBondBetweenAtoms(
                        new_bond_begin_idx, new_bond_end_idx
                    )
                    if bond.GetBondType() != new_bond.GetBondType():
                        sym_bond_entry = tuple(
                            sorted([bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()])
                        )
                        if sym_bond_entry not in ambiguous_bonds:
                            ambiguous_bonds.append(
                                tuple(
                                    sorted(
                                        [bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()]
                                    )
                                )
                            )
        if not ambiguous_bonds:
            ambiguous_atoms = (
                []
            )  # if no ambiguous bonds, there cannot be physically valid sets of ambiguous atoms
            # this is because that would imply that two different simple/graph connectivites can give different
            # formal charges, which is not supported in this implementation and likely not possible outside of
            # exotic transition metal groups
        return ambiguous_atoms, ambiguous_bonds
   
    def _substructures_disagree(smarts1, smarts2):
        ref1 = Chem.MolFromSmarts(smarts1)
        ref2 = Chem.MolFromSmarts(smarts2)

        Chem.SanitizeMol(
            ref1,
            Chem.SANITIZE_NONE,
        )
        Chem.SanitizeMol(
            ref2,
            Chem.SANITIZE_NONE,
        )

        Chem.SetAromaticity(ref1, Chem.AromaticityModel.AROMATICITY_MDL)
        Chem.SetAromaticity(ref2, Chem.AromaticityModel.AROMATICITY_MDL)

        # fuzzy1, neighbor_idxs = _fuzzy_query(ref1) 
        # sym_atoms, sym_bonds = _get_symmetrical_groups(fuzzy1, ref1)

        fuzzy2, neighbor_idxs = _fuzzy_query(ref2) 
        sym_atoms, sym_bonds = _get_symmetrical_groups(fuzzy2, ref2)

        # sym_atoms = list(set(sym_atoms1) & set(sym_atoms2))
        # sym_bonds = list(set(sym_bonds1) & set(sym_bonds2))

        for full_match in ref1.GetSubstructMatches(fuzzy2, maxMatches=0):

            for atom_i, j in zip(ref2.GetAtoms(), full_match):
                if atom_i.GetAtomicNum() == 0:  # ignore neighboring atoms
                    continue
                atom_j = ref1.GetAtomWithIdx(j)
                # error checking for overlapping substructures with priority. Enforce that no ambiguous
                # chemical assignments are made.
                if (
                    atom_i.GetFormalCharge() != atom_j.GetFormalCharge()
                    and atom_i.GetIdx() not in sym_atoms
                ):
                    return True

            for b in ref2.GetBonds():
                ref_bond_ids = tuple(
                    sorted([b.GetBeginAtomIdx(), b.GetEndAtomIdx()])
                )
                x = full_match[b.GetBeginAtomIdx()]
                y = full_match[b.GetEndAtomIdx()]
                b2 = ref1.GetBondBetweenAtoms(x, y)

                if (
                    b.GetBondType() != b2.GetBondType()
                    and ref_bond_ids not in sym_bonds
                ):
                    return True
        return False

    # check for disagreeing substructures in both directions
    if _substructures_disagree(s1, s2):
        return True
    elif _substructures_disagree(s2, s1):
        return True
    else:
        return False

def smarts_agree(smarts_input, json_dict):
    if len(json_dict) == 0:
        return True
    for monomer, smarts_list in json_dict.items():
        for smarts in smarts_list:
            if substructures_disagree(smarts, smarts_input):
                print(f"new query disagrees with {monomer}:\n{smarts}\n{smarts_input}")
                return False
    return True

GLOBAL_TOOLKIT_REGISTRY.deregister_toolkit(OpenEyeToolkitWrapper)

current_dir = Path(__file__).parent.resolve()
json_dir = current_dir / Path("json_files")
json_dir.mkdir(parents=False, exist_ok=True)
os.chdir(current_dir)

excluded_files = ["7wcc", "polyphenyleneI", "polyphenyleneII", "PolyphenyleneIII"] # files known to cause disagreements
# 7wcc includeds a FE+3 atom which disagrees with other forms or iron. Ions may need to be manually specified to load these files
# polyphenylene contains aromatic rings whose neighbor bonds disagree with the SN6 residue of 144d. 
json_dict = defaultdict(list)
for json_file in Path("json_files").glob("**/*.json"):
    if json_file.stem in excluded_files:
        continue
    print(json_file.stem)
    with open(json_file, "r") as file:
        data = json.load(file)
        monomers = data["monomers"]
        for mon, smarts_list in monomers.items():
            assert isinstance(smarts_list, list)
            for smarts in smarts_list:
                if smarts_agree(smarts, json_dict):
                    json_dict[json_file.stem].append(smarts)
                else:
                    print(f"\t^^disagreement for {json_file.stem}")
with open("combined_substruct_dict.json", "w") as file:
    json.dump(json_dict, file, indent=4)