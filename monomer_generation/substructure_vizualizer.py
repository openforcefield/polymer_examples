from rdkit import Chem
from pathlib import Path
from copy import deepcopy
from openff.toolkit.topology.molecule import Molecule
from openff.toolkit.topology.topology import Topology
# from openff.toolkit.utils.toolkits import RDKitToolkitWrapper, OpenEyeToolkitWrapper
from rdkit import Chem
from collections import defaultdict, OrderedDict
import py3Dmol
import ipywidgets as widgets
from IPython.display import display, clear_output
import time
import os
import tempfile
import json
import itertools

from openmm import unit as openmm_unit
from openmm.app import PDBFile
import numpy as np
from substructure_generator import SubstructureGenerator, Monomer
import matplotlib
from matplotlib import cm

import networkx as nx
from networkx.algorithms import isomorphism

class ChemistryEngine:
    def __init__(self, file):
        self.file = file
        self.substructure_generator = SubstructureGenerator()
        # make sure file exists
        path = Path(file)
        if not path.exists():
            print(f"file path {str(path)} does not exist, returning from class init")
            return
        suffix = path.suffix.lower()
        if 'pdb' in suffix:
            rdmol = Chem.rdmolfiles.MolFromPDBFile(self.file, removeHs=False, sanitize=False)
        elif 'sdf' in suffix:
            suppl = Chem.rdmolfiles.SDMolSupplier(self.file, removeHs=False, sanitize=False)
            rdmol = suppl[0]
        else:
            print("no valid file formats: only accepts pdb and sdf files")
        for atom in rdmol.GetAtoms():
            atom.SetAtomMapNum(atom.GetIdx())
            atom.SetProp("atomLabel", atom.GetSymbol() + f"{atom.GetAtomMapNum()}")
        self.full_molecule = rdmol
        self.n_atoms = rdmol.GetNumAtoms()
    
    def get_sdf_block(self, chemical_info={}):
        # additional_specs: dictionary of chemistry specs indexed by atom map number
        # example to get an sdf with added double bond:
        # block = self.get_sdf_block({'double': (21,23)})
        additional_specs = deepcopy(chemical_info)
        if "wildtypes" in additional_specs.keys():
            additional_specs.pop("wildtypes")
        mol_copy = self.assign_additional_specs(additional_specs)
    
        with tempfile.TemporaryDirectory() as tmpdir:
            prev_dir = os.getcwd()
            os.chdir(os.path.abspath(tmpdir))
            writer = Chem.rdmolfiles.SDWriter('molecule.sdf')
            writer.write(mol_copy)

            with open("molecule.sdf", "r") as file:
                sdf_block = file.read()
            os.chdir(prev_dir)
        return sdf_block

    def assign_additional_specs(self, additional_specs):
        mol_copy = deepcopy(self.full_molecule)
        # modify the molecule:
        doubles = additional_specs.get('double', [])
        sorted_doubles = [tuple(sorted(t)) for t in doubles]

        triples = additional_specs.get('triple', [])
        sorted_triples = [tuple(sorted(t)) for t in triples]

        for bond in mol_copy.GetBonds():
            bond_ids = tuple(sorted([bond.GetBeginAtom().GetAtomMapNum(), bond.GetEndAtom().GetAtomMapNum()]))
            if bond_ids in sorted_doubles:
                bond.SetBondType(Chem.rdchem.BondType.DOUBLE)
            elif bond_ids in sorted_triples:
                bond.SetBondType(Chem.rdchem.BondType.TRIPLE)

        wildtypes = additional_specs.get('wildtypes', [])
        for atom in mol_copy.GetAtoms():
            if atom.GetIdx() in wildtypes:
                atom.SetAtomicNum(0)

        return mol_copy

    def get_wildtypes(self, ids):
        wildtypes = set()
        for atom in self.full_molecule.GetAtoms():
            if atom.GetIdx() in ids:
                # check all neighbors
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetIdx() not in ids:
                        wildtypes.add(neighbor.GetIdx())
        return list(wildtypes)
    
    def get_chemical_info(self, rdmol):
        # modify the molecule:
        chemical_info = defaultdict(list)

        for bond in rdmol.GetBonds():
            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                # bond_ids = tuple(sorted([bond.GetBeginAtom().GetAtomMapNum(), bond.GetEndAtom().GetAtomMapNum()]))
                bond_ids = tuple(sorted([bond.GetBeginAtom().GetIdx(), bond.GetEndAtom().GetIdx()]))
                chemical_info["double"].append(bond_ids)
            elif bond.GetBondType() == Chem.rdchem.BondType.TRIPLE:
                bond_ids = tuple(sorted([bond.GetBeginAtom().GetIdx(), bond.GetEndAtom().GetIdx()]))
                chemical_info["triple"].append(bond_ids)

        for atom in rdmol.GetAtoms():
            if atom.GetAtomicNum() == 0:
                chemical_info['wildtypes'].append(atom.GetIdx())

        return chemical_info
    
    def update_chemical_info_ids(self, isomorphism, chemical_info):
        new_chemical_info = defaultdict(list)
        for key, values in chemical_info.items():
            for value in values:
                if isinstance(value, int):
                    new_chemical_info[key].append(isomorphism[value])
                elif isinstance(value, tuple):
                    i, j = value
                    new_i = isomorphism[i]
                    new_j = isomorphism[j]
                    new_chemical_info[key].append(tuple([new_i, new_j]))
        return new_chemical_info

    def get_smarts_from_ids(self, ids, additional_specs):
        mol_copy = self.assign_additional_specs(additional_specs)
        # create atommap numbers and a mapping dict to be used with this monomer
        map_num = 1
        ids_to_map_num = {}
        for atom in mol_copy.GetAtoms():
            if atom.GetIdx() in ids and atom.GetIdx() not in additional_specs.get('wildtypes', []):
                atom.SetAtomMapNum(map_num)
                ids_to_map_num[atom.GetIdx()] = map_num
                map_num += 1
            else:
                atom.SetAtomMapNum(0)
        ids += additional_specs.get('wildtypes', []) # update ids to include any caps
        monomer_smarts = Chem.MolFragmentToSmarts(mol_copy, atomsToUse = ids)
        # must read and write again to update smarts query to print "*" instead of "#0" -> TODO: why tf does this happen?
        rdmol = Chem.MolFromSmarts(monomer_smarts)
        for atom in rdmol.GetAtoms():
            if atom.GetAtomicNum() == 0:
                atom.SetQuery(Chem.AtomFromSmarts("[*]"))
        monomer_smarts = Chem.MolToSmarts(rdmol)
        # finally, return a new "additional_specs/chemical_info" dict
        chemical_info = self.get_chemical_info(rdmol)
        return monomer_smarts, rdmol, chemical_info

    def add_monomer_from_ids(self, name, ids, additional_specs):
        monomer_smarts, _, _ = self.get_smarts_from_ids(ids, additional_specs)
        self.substructure_generator.add_monomer(name, monomer_smarts)

    def test_polymer_load(self, substructure_lib, verbose=True):
        # executes from_pdb using the given biopolymer substructure library and records how many
        # atoms and bonds have chemical information assigned.
        # returns: original molecule ids that are assigned and have assigned bonds 
        assigned_atoms = set()
        unassigned_atoms = set()
        assigned_bonds = set()
        unassigned_bonds = set()

        topology, isomorphism_summary, error = Topology.from_pdb_and_monomer_info(self.file, substructure_lib, strict=False)
        chemical_info = {"double": [], "triple": []}
        for atom in topology.atoms:
            if atom.metadata['already_matched']:
                assigned_atoms.add(atom.molecule_atom_index)
            else:
                unassigned_atoms.add(atom.molecule_atom_index)
        for bond in topology.bonds:
            # check for assigned bonds 
            if bond.bond_order in [1,1.5,2,3]:
                assigned_bonds.add((bond.atom1_index, bond.atom2_index))
                if bond.bond_order == 2:
                    chemical_info["double"].append((bond.atom1_index, bond.atom2_index))
                if bond.bond_order == 3:
                    chemical_info["triple"].append((bond.atom1_index, bond.atom2_index))
            else:
                unassigned_bonds.add((bond.atom1_index, bond.atom2_index))
                
        # print info on the polymer loading
        if verbose:
            print(f"number of atoms assigned: {len(assigned_atoms)}")
            print(f"number of bonds assigned: {len(assigned_bonds)}")
            print(f"number of atoms not assigned: {len(unassigned_atoms)}")
            print(f"number of bonds not assigned: {len(unassigned_bonds)}")
        return assigned_atoms, assigned_bonds, unassigned_atoms, unassigned_bonds, chemical_info, isomorphism_summary

    def _pdb_to_networkx(self, pdb_file):
        """
        Construct an OpenFF Topology object from an OpenMM Topology object.

        Parameters
        ----------
        substructure_library : dict{str:list[str, list[str]]}
            A dictionary of substructures. substructure_library[aa_name] = list[tagged SMARTS, list[atom_names]]
        openmm_topology : openmm.app.Topology
            An OpenMM Topology object

        Returns
        -------
        omm_topology_G : networkx graph
            A networkX graph representation of the openmm topology with chemical information added from the
            substructure dictionary. Atoms are nodes and bonds are edges.
            Nodes (atoms) have attributes for `atomic_number` (int) and `formal_charge` (int).
            Edges (bonds) have attributes for `bond_order` (Chem.rdchem.BondType).
            Any edges that are not assgined a bond order will have the value Chem.rdchem.BondType.UNSPECIFIED
            and should be considered an error.
        """
        import networkx as nx
        from openmm.app import PDBFile

        pdb = PDBFile(pdb_file)
        openmm_topology = pdb.topology

        omm_topology_G = nx.Graph()
        for atom in openmm_topology.atoms():
            omm_topology_G.add_node(
                atom.index,
                atomic_number=atom.element.atomic_number,
                formal_charge=0.0,
                atom_name=atom.name,
                residue_name=atom.residue.name,
                residue_number=atom.residue.index,
            )

        n_hydrogens = [0] * openmm_topology.getNumAtoms()
        for bond in openmm_topology.bonds():
            omm_topology_G.add_edge(
                bond.atom1.index,
                bond.atom2.index,
                bond_order=Chem.rdchem.BondType.UNSPECIFIED,  # bond.order
            )
            # omm_topology_G.add_edge(
            #     bond.atom1.index,
            #     bond.atom2.index,
            #     bond_order=1,  # quick fix for now
            # )
            # Assign sequential negative numbers as atomic numbers for hydrogens attached to the same heavy atom.
            # We do the same to the substructure templates that are used for matching. This saves runtime because
            # it removes redundant self-symmetric matches.
            if bond.atom1.element.atomic_number == 1:
                h_index = bond.atom1.index
                heavy_atom_index = bond.atom2.index
                n_hydrogens[heavy_atom_index] += 1
                omm_topology_G.nodes[h_index]["atomic_number"] = (
                    -1 * n_hydrogens[heavy_atom_index]
                )
            if bond.atom2.element.atomic_number == 1:
                h_index = bond.atom2.index
                heavy_atom_index = bond.atom1.index
                n_hydrogens[heavy_atom_index] += 1
                omm_topology_G.nodes[h_index]["atomic_number"] = (
                    -1 * n_hydrogens[heavy_atom_index]
                )

        # Try matching this substructure to the whole molecule graph
        # node_match = isomorphism.categorical_node_match(
        #     ["atomic_number", "already_matched"], [-100, False]
        # )
        return omm_topology_G

    def _rdmol_to_networkx(self, rdmol):
        _bondtypes = {
            # 0: Chem.BondType.AROMATIC,
            Chem.BondType.SINGLE: 1,
            Chem.BondType.AROMATIC: 1.5,
            Chem.BondType.DOUBLE: 2,
            Chem.BondType.TRIPLE: 3,
            Chem.BondType.QUADRUPLE: 4,
            Chem.BondType.QUINTUPLE: 5,
            Chem.BondType.HEXTUPLE: 6,
        }
        rdmol_G = nx.Graph()
        n_hydrogens = [0] * rdmol.GetNumAtoms()
        for atom in rdmol.GetAtoms():
            atomic_number = atom.GetAtomicNum()
            # Assign sequential negative numbers as atomic numbers for hydrogens attached to the same heavy atom.
            # We do the same to hydrogens in the protein graph. This makes it so we
            # don't have to deal with redundant self-symmetric matches.
            if atomic_number == 1:
                heavy_atom_idx = atom.GetNeighbors()[0].GetIdx()
                n_hydrogens[heavy_atom_idx] += 1
                atomic_number = -1 * n_hydrogens[heavy_atom_idx]
            
            rdmol_G.add_node(
                atom.GetIdx(),
                atomic_number=atomic_number,
                formal_charge=atom.GetFormalCharge(),
                map_num=atom.GetAtomMapNum()
            )
            # These substructures (and only these substructures) should be able to overlap previous matches.
            # They handle bonds between substructures.
        for bond in rdmol.GetBonds():
            bond_type = bond.GetBondType()

            # All bonds in the graph should have been explicitly assigned by this point.
            if bond_type == Chem.rdchem.BondType.UNSPECIFIED:
                raise Exception
                # bond_type = Chem.rdchem.BondType.SINGLE
                # bond_type = Chem.rdchem.BondType.AROMATIC
                # bond_type = Chem.rdchem.BondType.ONEANDAHALF
            rdmol_G.add_edge(
                bond.GetBeginAtomIdx(),
                bond.GetEndAtomIdx(),
                bond_order=_bondtypes[bond_type],
            )
        return rdmol_G

    def _get_isomorphisms(self, query, structure):
        # returns an isomorphism map using networkx from query to structure where 
        # both query and structure are rdkit molecules 

        def node_match(data1, data2):
            if data1.get("atomic_number", -100) == 0 or data2.get("atomic_number", -100) == 0:
                return True
            elif data1.get("atomic_number", -100) == data2.get("atomic_number", -100):
                if data1.get("atom_map", 0) == data2.get("atom_map", 0):
                    # return true if atomic numbers match on non_captured atoms
                    return True
                else:
                    # else, captured atoms must have maching values for "already_matched"
                    if data1.get("already_matched", False) == data2.get("already_matched", False):
                        return True
                    else:
                        return False
            else:
                return False

        if isinstance(structure, str):
            if ".pdb" in structure or ".PDB" in structure:
                rdmol_G = self._rdmol_to_networkx(query)
                omm_topology_G = self._pdb_to_networkx(structure)
                GM = isomorphism.GraphMatcher(
                    omm_topology_G, rdmol_G, node_match=node_match
                )
                return GM.subgraph_is_isomorphic(), GM.subgraph_isomorphisms_iter()
            else:
                return -1, -1
        elif isinstance(structure, Chem.rdchem.Mol):
            rdmol_G = self._rdmol_to_networkx(query)
            structure_G = self._rdmol_to_networkx(structure)
            GM = isomorphism.GraphMatcher(
                structure_G, rdmol_G, node_match=node_match
            )
            return GM.subgraph_is_isomorphic(), GM.subgraph_isomorphisms_iter()
        else:
            return -1, -1

    def find_terminal_groups(self, ids, additional_specs):
        monomer_smarts, monomer, chemical_info = self.get_smarts_from_ids(ids, additional_specs)
        open_atoms_query = chemical_info["wildtypes"]
        n_atoms = monomer.GetNumAtoms() - len(open_atoms_query)
        is_isomorphic, isomorphisms = self._get_isomorphisms(monomer, self.full_molecule)
        mapped_atoms_set = set()
        terminal_groups = []
        
        if is_isomorphic:
            isomorphisms_list = list(isomorphisms)
            for iso in isomorphisms_list:
                keys = set(iso.keys())
                for struct_id, query_id in iso.items():
                    if query_id in open_atoms_query:
                        keys.remove(struct_id)
                mapped_atoms_set = mapped_atoms_set | keys

            for iso, id in zip(isomorphisms_list, range(0, len(isomorphisms_list))):
                iso_chemical_info = deepcopy(chemical_info)
                substructure_ids = set(iso.keys()) # start with the initial isomorphism
                original_ids = set(iso.keys())
                open_atoms = []
                for struct_id, query_id in iso.items():
                    if query_id in open_atoms_query:
                        open_atoms.append(struct_id)
                terminal_group = False
                for open_atom in open_atoms:
                    queue = [open_atom]
                    visited_atoms = set()
                    searching = True
                    isomorphism_adjacent = False
                    end_reached = False
                    while searching:
                        current_id = queue[0]
                        current_atom = self.full_molecule.GetAtomWithIdx(current_id)
                        for neighbor in current_atom.GetNeighbors():
                            idx = neighbor.GetIdx()
                            if idx not in visited_atoms and idx not in original_ids:
                                queue.append(idx)
                        queue.pop(0)
                        visited_atoms.add(current_id)
                        # now test for if the search can stop
                        # condition 1: the end of the graph is more than 20 atoms away:
                        if len(visited_atoms) > n_atoms:
                            searching = False
                        # condition 2: the end of the graph is reached
                        elif len(queue) == 0:
                            searching = False
                            terminal_group = True
                            end_reached = True
                    if end_reached:
                        substructure_ids = substructure_ids | visited_atoms
                        iso_chemical_info["wildtypes"].remove(iso[open_atom])
                if terminal_group:
                    reverse_iso = dict([(j,i) for i,j in iso.items()])
                    new_chemical_info = self.update_chemical_info_ids(reverse_iso, iso_chemical_info)
                    wildtypes = new_chemical_info.get("wildtypes", [])
                    substructure_ids = substructure_ids - set(wildtypes)
                    terminal_groups.append((list(substructure_ids), new_chemical_info))
        else:
            return []
        return terminal_groups

class PolymerVisualizer3D:
    instances = {}
    def __init__(self, chemistry_engine):
        self.instance_id = int(time.time() * 1000.0)
        self.__class__.instances[self.instance_id] = self
        self.chemistry_engine = chemistry_engine
        sdf_block = chemistry_engine.get_sdf_block()
        self.width = 800
        self.height = 500
        self.selections = defaultdict(Selection)
        self.highlights = set()
        self.view = self._view_from_sdf_block(sdf_block)
        self.view.setStyle({"model": -1}, {"stick": {}})
        # buttons
        next_button = widgets.Button(description="Next")
        prev_button = widgets.Button(description="Previous")
        finish_button = widgets.Button(description="Finish Monomer")
        double_bonds_button = widgets.Button(description="Make Double bonds")
        triple_bonds_button = widgets.Button(description="Make Triple bonds")
        formal_charge_button = widgets.Button(description="assign formal charges")
        formal_charge_menu = widgets.Dropdown(
                                                options=['-3','-2','-1','0','1','2','3'],
                                                value='0',
                                                disabled=False,
                                                description='',
                                                layout = widgets.Layout(width='50px')
                                            )
        run_button = widgets.Button(description="Run")
        print_button = widgets.Button(description="Print Selection")
        color_code_button = widgets.Button(description="Inspect", disabled=True)
        remove_highlights_button = widgets.Button(description="Remove Red Highlights", disabled=True)
        header_tags = widgets.ToggleButtons(
                                        options=['Test Load', 'Edit Monomers'],
                                        description='',
                                        disabled=False,
                                        button_style='', # 'success', 'info', 'warning', 'danger' or ''
                                        tooltips=[''],
                                        value='Test Load',
                                        layout = widgets.Layout(width=f'{self.width-1}px'),
                                        style={"button_width": f'{self.width/2 - 10}px'}
                                    )
        editing_tags = widgets.ToggleButtons(
                                        options=['Select Atoms', 'Assign Chemical Info', 'Inspect Caps'],
                                        description='',
                                        disabled=True,
                                        button_style='', # 'success', 'info', 'warning', 'danger' or ''
                                        tooltips=[''],
                                        value='Select Atoms',
                                        layout = widgets.Layout(width='auto')
                                    )
        terminal_group_tags = widgets.ToggleButtons(
                                        options=["No Terminal Groups Found"],
                                        description='',
                                        disabled=True,
                                        button_style='', # 'success', 'info', 'warning', 'danger' or ''
                                        tooltips=[''],
                                        value='No Terminal Groups Found',
                                        layout = widgets.Layout(width='auto')
        )
        name_input_box = widgets.Text(
                                        value='',
                                        placeholder='Monomer Name',
                                        description='',
                                        disabled=False,
                                        layout = widgets.Layout(width='150px')
                                    )
        terminal_group_delete_button = widgets.Button(description="Del Terminal Group")
        selected_atom_label = widgets.Label("")
        iso_inspect_dropdown = widgets.Dropdown(
                                                    options=['Fitted Isomorphisms', 'Unmapped Atoms'],
                                                    value='Fitted Isomorphisms',
                                                    description='',
                                                    disabled=True,
                                                )
        iso_inspect_dropdown.style.description_width = 'auto'                                        
        color_box = widgets.ColorPicker(
                                            concise=True,
                                            description='',
                                            value='#FFFFFF',
                                            disabled=True
                                        )
        valid_check = widgets.Valid(
                                        value=False,
                                        description=''
                                    )
        valid_check.readout = "Click \"Run\""
        
        next_button.on_click(self._function_next_button)
        prev_button.on_click(self._function_prev_button)
        finish_button.on_click(self._function_finish_button)
        double_bonds_button.on_click(self._function_double_bonds_button)
        triple_bonds_button.on_click(self._function_triple_bonds_button)
        formal_charge_button.on_click(self._function_formal_charge_button)
        run_button.on_click(self._function_run_button)
        print_button.on_click(self._function_print_button)
        color_code_button.on_click(self._function_color_code_button)
        remove_highlights_button.on_click(self._function_remove_highlights_button)
        header_tags.observe(self._function_header_tags, "value")
        terminal_group_tags.observe(self._function_terminal_group_tags, "value")
        terminal_group_delete_button.on_click(self._function_terminal_group_delete_button)
        iso_inspect_dropdown.observe(self._function_iso_inspect_dropdown, "value")

        # any references to any widgets should only be made to this dict 
        self.widgets = {
                        "next_button": next_button,
                        "prev_button": prev_button, 
                        "finish_button": finish_button,
                        "double_bonds_button": double_bonds_button, 
                        "triple_bonds_button": triple_bonds_button, 
                        "formal_charge_button": formal_charge_button,
                        "formal_charge_menu": formal_charge_menu,
                        "run_button": run_button, 
                        "print_button": print_button,
                        "color_code_button": color_code_button, 
                        "remove_highlights_button": remove_highlights_button,
                        "header_tags": header_tags, 
                        "editing_tags": editing_tags,
                        "terminal_group_tags": terminal_group_tags,
                        "name_input_box": name_input_box,
                        "terminal_group_delete_button": terminal_group_delete_button,
                        "selected_atom_label": selected_atom_label,
                        "iso_inspect_dropdown": iso_inspect_dropdown,
                        "color_box": color_box,
                        "valid_check": valid_check}
        
        self.single_selected_atom = -1       
        self.menu_mode = "test_load-default"
        self.valid_menu_modes = ['test_load-default', 'edit-select_monomer', 'edit-assign_chemistry', 'edit-inspect_caps']
        self.click_mode = "select_monomer"
        self.valid_click_modes = ['select_monomer', 'do_nothing', 'select_double', 'select_triple', 'select_charged']

        # after isomorphisms are generated:
        self.isomorphisms = []
        self.monomer_colors = dict()

    # the only function the user should ever have to call themselves 
    def show(self):
        self._reset_view()
        return

    def __repr__(self):
        self.show()
        return f"file: {self.chemistry_engine.file}"

    def _add_highlights(self, ids):
        for atom in ids:
            self.view.addSphere({"center":{"serial": atom},
                                "radius": 0.75,
                                "color": 'red',
                                "opacity": "0.3"})
        self.highlights = ids
    
    def _view_from_sdf_block(self, block):
        view = py3Dmol.view(width=self.width, height=self.height)
        view.addModel(block, "sdf", {"keepH": True})
        view.center({"model": -1})
        return view

    #MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
    #-------------------------------Clickables-------------------------------
    #WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
    def _set_clickable(self, new_mode):
        if new_mode not in self.valid_click_modes:
            print("internal: new_mode is not valid")
            return
        self.click_mode = new_mode
        if new_mode == "select_monomer":
            selection = {"model": -1}
            code = '''function(atom,viewer,event,container) {
                    if(!atom.style['clicksphere']) {
                            viewer.setStyle({"serial": atom.serial}, {"stick": {"color": 0xFF0000}, "clicksphere": {"radius": 0.25}})
                            
                            var serial = atom.serial;
                            Jupyter.notebook.kernel.execute("PolymerVisualizer3D.instances[%d]._store_selected_atom(" + serial + ")");
                    }
                    else {
                            viewer.setStyle({"serial": atom.serial}, {"stick": {"colorscheme": "default"}})
                        
                            var serial = atom.serial;
                            Jupyter.notebook.kernel.execute("PolymerVisualizer3D.instances[%d]._remove_selected_atom(" + serial + ")");
                    }
                    viewer.render();}''' % (self.instance_id, self.instance_id)
        elif new_mode == 'do_nothing':
            selection = {"model": -1}
            code = '''function(atom,viewer,event,container) {
                        var serial = atom.serial;
                        Jupyter.notebook.kernel.execute("PolymerVisualizer3D.instances[%d]._store_selected_atom(" + serial + ")");
                    }''' % (self.instance_id)
        elif new_mode in ["select_double", "select_triple", "select_captured"]:
            selection_name = self._get_selection()
            selection = {"serial": self.selections[selection_name].selected_atoms['select_monomer'] + self.selections[selection_name].chemical_info.get("wildtypes", [])}
            # selection = {"model": -1}
            code = '''function(atom,viewer,event,container) {
                        if(!atom.style['clicksphere']) {
                                viewer.setStyle({"serial": atom.serial}, {"stick": {"color": 0xFFFF00}, "clicksphere": {"radius": 0.25}})
            
                                var serial = atom.serial;
                                Jupyter.notebook.kernel.execute("PolymerVisualizer3D.instances[%d]._store_selected_atom(" + serial + ")");
                        }
                        else {
                                viewer.setStyle({"serial": atom.serial}, {"stick": {"colorscheme": "default"}})
     
                                var serial = atom.serial;
                                Jupyter.notebook.kernel.execute("PolymerVisualizer3D.instances[%d]._remove_selected_atom(" + serial + ")");
                        }
                        viewer.render();}''' % (self.instance_id, self.instance_id)
        elif new_mode == "select_charged":
            selection_name = self._get_selection()
            selection = {"serial": self.selections[selection_name].selected_atoms['select_monomer']}
            return #not yet implemented
        else:
            print("internal error in _set_clickable")
            return

        self.view.setClickable(selection,True,code)
    
    def _store_selected_atom(self, serial):

        # print(f"in the store with serial number: {serial}")
        if self.click_mode == 'do_nothing':
            self.single_selected_atom = serial
            self.widgets["selected_atom_label"].value = self.get_atom_description(serial)
            return
        selection_name = self._get_selection()
        data = self.selections[selection_name].selected_atoms[self.click_mode]
        if serial not in data:
            self.selections[selection_name].selected_atoms[self.click_mode].append(serial)
        else:
            print("error: serial number already in self.selections[\"main\"].selected_atoms")

    def _remove_selected_atom(self, serial):
        selection_name = self._get_selection()
        data = self.selections[selection_name].selected_atoms[self.click_mode]
        if serial in data:
            self.selections[selection_name].selected_atoms[self.click_mode].remove(serial)
        else:
            print("error: serial not in self.selections[\"main\"].selected_atoms")
    
    def get_atom_description(self, serial):
        monomer_name = ""
        for name, ids, selected in self.isomorphisms:
            if serial in ids:
                monomer_name = name
                if selected:
                    break

        if monomer_name != "":
            return f"{monomer_name} atom with serial: {serial}"
        else:
            return f"unmatched atom with serial: {serial}"
    #MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
    #-------------------------------Hoverables-------------------------------
    #WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
    def _set_hoverable(self, mode):
        if mode not in self.valid_click_modes:
            print("internal: new_mode is not valid")
            return
        func1 = '''
            function(atom,viewer,event,container) {
                    if(!atom.label) {
                        atom.label = viewer.addLabel(atom.elem+":"+ atom.serial,{position: atom, backgroundColor: 'mintcream', fontColor:'black', backgroundOpacity: "0.3"});
                    }
                }
            '''
        func2 = '''
                function(atom,viewer,event,container) {
                        if(atom.label) {
                            viewer.removeLabel(atom.label);
                            delete atom.label;
                        }
                    }
            '''
        # if mode != 'select_monomer':
        #     selection = {'serial': self.selections["main"].selected_atoms['select_monomer'], 'invert': True}
        # else:
        #     selection = {'model': -1}
        selection = {'model': -1}
        self.view.setHoverable(selection, True, func1, func2)
        self.view.setHoverDuration(100)

    #MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
    #------------------------------View Actions------------------------------
    #WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
    def _reload_view(self, center_to_monomer=False):
        # called anytime the view is changed
        clear_output(wait=False)
        display(self.view, self.buttons)
    
    def _reset_view(self):
        # removes shapes and resets all clickables/views
        self.view.setStyle({"model": -1}, {"stick": {"colorscheme": "default"}})
        self.view.removeAllLabels()
        self.selections = defaultdict(Selection)
        self.selections["main"] = Selection()
        self._set_clickable('do_nothing')
        self._set_hoverable('do_nothing')
        self._load_buttons('test_load-default')
        self._reload_view()

    def _set_style_monomer_view(self, cap=""):
        selection_name = self._get_selection()
        if cap == "":
            body_ids = self.selections[selection_name].selected_atoms['select_monomer']
            cap_ids = self.selections[selection_name].chemical_info['wildtypes']
        else:
            body_ids = self.terminal_groups[cap].selected_atoms['select_monomer']
            cap_ids = self.terminal_groups[cap].chemical_info['wildtypes']
        sdf_block = self.chemistry_engine.get_sdf_block(self.selections[selection_name].chemical_info)
        self.highlights = set()
        self.view = self._view_from_sdf_block(sdf_block)
        self.view.removeAllShapes()
        self.view.center({"serial": body_ids + cap_ids})
        self.view.setStyle({"serial": body_ids}, {"stick": {"colorscheme": "default"}})
        self.view.setStyle({"serial": cap_ids}, {"stick": {"color": "0xFFFF00"}})
        self.view.setStyle({"serial": body_ids + cap_ids, "invert": True}, {"line": {}})
    #MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
    #------------------------------Button Utils------------------------------
    #WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
    def _load_buttons(self, new_menu_mode=""):
        spacer = widgets.Label("")
        header_tag_layout = widgets.Layout(display='flex',
                                        align_items='center',
                                        justify_content='center',
                                        width=f'{self.width}px')
        editing_tag_layout = widgets.Layout(display='flex',
                                        align_items='center',
                                        justify_content='center',
                                        width=f'{self.width}px')
        header_tags = widgets.HBox(children=[self.widgets["header_tags"]], layout=header_tag_layout)
        editing_tags = widgets.HBox(children=[self.widgets["editing_tags"]], layout=editing_tag_layout)
        if new_menu_mode != "":
            self.menu_mode = new_menu_mode
        if self.menu_mode == 'test_load-default':
            self.widgets['header_tags'].disabled = False
            self.widgets['color_code_button'].disabled = True
            self.widgets['valid_check'].value = False
            self.widgets['valid_check'].readout = "Click \"Run\""
            run_button_box = widgets.HBox((self.widgets['run_button'], self.widgets['valid_check']))
            color_codes_box = widgets.HBox((self.widgets['color_code_button'], self.widgets['iso_inspect_dropdown'], self.widgets['color_box']))
            self.buttons = widgets.VBox((header_tags, run_button_box, self.widgets["print_button"], color_codes_box, self.widgets["selected_atom_label"]))
        
        elif self.menu_mode == 'edit-select_monomer':
            self.widgets['editing_tags'].value = "Select Atoms"
            self.widgets['header_tags'].disabled = False
            if len(self.highlights) > 0:
                self.widgets['remove_highlights_button'].disabled = False
            else:
                self.widgets['remove_highlights_button'].disabled = True
            description = widgets.Label("Select atoms to build a monomer unit and press \"Next\". Be sure to select one contiguous set of atoms.")

            self.buttons = widgets.VBox((header_tags, editing_tags, description, self.widgets["remove_highlights_button"], spacer, self.widgets["next_button"]))

        elif self.menu_mode == 'edit-assign_chemistry':
            self.widgets['editing_tags'].value = "Assign Chemical Info"
            self.widgets['header_tags'].disabled = True
            description = widgets.Label("Assign chemistry to the monomer. Click a button once to initiate editing mode.")
            c1 = widgets.VBox((self.widgets['double_bonds_button'], self.widgets['triple_bonds_button']))
            c2 = widgets.VBox((self.widgets['formal_charge_button'], self.widgets['formal_charge_menu']))
            chemistry_buttons = widgets.HBox((c1,c2))
            prev_next_buttons = widgets.HBox((self.widgets["prev_button"], self.widgets["next_button"]))

            self.buttons = widgets.VBox((header_tags, editing_tags, description, chemistry_buttons, spacer, prev_next_buttons))

        elif self.menu_mode == 'edit-inspect_caps':
            self.widgets['editing_tags'].value = "Inspect Caps"
            self.widgets['header_tags'].disabled = True
            description = widgets.Label("Inspect all caps that could be found in the pdb (and assign chemical info if needed). All caps not found in the pdb must be added manually")
            c1 = widgets.VBox((self.widgets['double_bonds_button'], self.widgets['triple_bonds_button']))
            c2 = widgets.VBox((self.widgets['formal_charge_button'], self.widgets['formal_charge_menu']))
            chemistry_buttons = widgets.HBox((c1,c2))
            prev_next_buttons = widgets.HBox((self.widgets["prev_button"], self.widgets["name_input_box"], self.widgets["finish_button"]))
            # create buttons for any selections other than "main"
            button_values = []
            for key in self.selections.keys():
                if key != "main":
                    button_values.append(key)
            if len(button_values) > 0:
                self.widgets["terminal_group_tags"].options = button_values
                self.widgets["terminal_group_tags"].value = button_values[0]
                self.widgets["terminal_group_tags"].disabled = False
            else:
                self.widgets["terminal_group_tags"].options = ["No Terminal Groups Found"]
                self.widgets["terminal_group_tags"].value = "No Terminal Groups Found"
                self.widgets["terminal_group_tags"].disabled = True
            terminal_groups_row = widgets.HBox((self.widgets["terminal_group_tags"], self.widgets["terminal_group_delete_button"]))

            self.buttons = widgets.VBox((header_tags, editing_tags, description, terminal_groups_row, chemistry_buttons, spacer, prev_next_buttons))

            self._set_style_monomer_view()
        else:
            print("invalid new menu mode")

    def _enable_buttons(self, ids):
        for id in self.widgets.keys():
            if id in ids:
                self.widgets[id].disabled = False
    def _disable_buttons(self, ids):
        for id in self.widgets.keys():
            if id in ids:
                self.widgets[id].disabled = True

    def _get_selection(self):
        if self.menu_mode == 'edit-inspect_caps':
            selection = self.widgets['terminal_group_tags'].value
            if selection == "No Terminal Groups Found":
                print("no terminal groups to delete")
        else:
            selection = "main"
        return selection
    #MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
    #-------------------------Button Click Actions---------------------------
    #WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
    def _function_next_button(self, b):
        if self.widgets['editing_tags'].value == "Select Atoms":
            # load the selection
            self.selections["main"].chemical_info['wildtypes'] = self.chemistry_engine.get_wildtypes(self.selections["main"].selected_atoms['select_monomer'])
            self._set_style_monomer_view()
            self._load_buttons("edit-assign_chemistry")
            
        elif self.widgets['editing_tags'].value == "Assign Chemical Info":
            # find isomorphisms at the end of the molecule using chemistry_engine
            main_selected_ids = self.selections["main"].selected_atoms["select_monomer"]
            main_chemical_info = self.selections["main"].chemical_info
            terminal_group_isomorphisms = self.chemistry_engine.find_terminal_groups(main_selected_ids, main_chemical_info)
            terminal_group_name = "TERM"
            terminal_group_name_id = 1
            for selected_ids, chemical_info in terminal_group_isomorphisms:
                terminal_group_selection = Selection()
                terminal_group_selection.selected_atoms = {"select_monomer": selected_ids}
                terminal_group_selection.chemical_info = chemical_info
                self.selections[f"{terminal_group_name}{terminal_group_name_id}"] = terminal_group_selection
                terminal_group_name_id += 1
                pass
            self._load_buttons("edit-inspect_caps")
        else:
            print("next button error")
        self._set_clickable('do_nothing')
        self._set_hoverable('do_nothing')
        self._reload_view()

    def _function_prev_button(self, b):
        if self.widgets['editing_tags'].value == "Inspect Caps":
            main_selection = self.selections["main"]
            self.selections = defaultdict(Selection)
            self.selections["main"] = main_selection
            self._load_buttons("edit-assign_chemistry")
            
            self._set_style_monomer_view()
            self._set_clickable('do_nothing')
            self._set_hoverable('do_nothing')
        elif self.widgets['editing_tags'].value == "Assign Chemical Info":
            self._load_buttons("edit-select_monomer")
            self._set_clickable('select_monomer')
            self._set_hoverable('select_monomer')

            monomer_atoms = self.selections["main"].selected_atoms['select_monomer']
            self.selections["main"].selected_atoms = defaultdict(list)
            self.selections["main"].selected_atoms[self.click_mode] = monomer_atoms
            self.selections["main"].chemical_info = defaultdict(list)
            self.view.center({"serial": self.selections["main"].selected_atoms['select_monomer']})
            self.view.setStyle({"serial": self.selections["main"].selected_atoms['select_monomer']}, {"stick": {"color": 0xFF0000}, "clicksphere": {"radius": 0.25}})
            self.view.setStyle({"serial": self.selections["main"].selected_atoms['select_monomer'], "invert": True}, {"stick": {"colorscheme": "default"}})
        else:
            print("prev button error")
        
        self._reload_view()
        
    def _function_finish_button(self, b):
        if self.widgets['editing_tags'].value != "Inspect Caps":
            print("finish button error")

        # create a monomer entry for the SubstructureGenerator Class through ChemistryEngine
        name = self.widgets['name_input_box'].value
        for selection_name in self.selections.keys():
            atoms = deepcopy(self.selections[selection_name].selected_atoms['select_monomer'])
            if len(atoms) == 0:
                continue 
            info = deepcopy(self.selections[selection_name].chemical_info)
            if selection_name == "main":
                monomer_name = name
            else:
                monomer_name = name + "_" + selection_name
            self.chemistry_engine.add_monomer_from_ids(monomer_name, atoms, info)
        
        self._load_buttons("edit-select_monomer")
        self._set_clickable('select_monomer')
        self._set_hoverable('select_monomer')

        self.selections = defaultdict(Selection)
        self.selections["main"] = Selection()
        self.view.removeAllShapes()
        self.view.setStyle({"model": -1}, {"stick": {"colorscheme": "default"}})
 
        self._reload_view()
        
    def _function_double_bonds_button(self, b):
        temp_disable_widgets = ["triple_bonds_button", "formal_charge_button", "formal_charge_menu"]
    
        # check where to store info, either in self.selections["main"] or in one of self.terminal_groups
        selection = self._get_selection()
            
        # cycles between "Make Double bonds" and "Commit Double Bond"
        if b.description == "Make Double bonds":
            b.description = "Commit Double Bond"
            b.button_style = "warning"
            self._set_clickable('select_double')
            self._set_hoverable('select_double')
            self.view.center({"serial": self.selections[selection].selected_atoms['select_monomer']})
            self._disable_buttons(temp_disable_widgets)
            self._reload_view()
        elif b.description == "Commit Double Bond":
            b.description = "Make Double bonds"
            b.button_style = ""
            
            # load the new double bond into the view with the chemistry engine + sdf
            # check that only two atoms were selected
            bond = self.selections[selection].selected_atoms['select_double']
            self.selections[selection].selected_atoms['select_double'] = []
            if len(bond) == 2:
                new_bond = tuple(bond) # sdf format indexed from 0
                bonds = self.selections[selection].chemical_info['double'] # does default dict handle this? TODO
                bonds.append(new_bond)
                self.selections[selection].chemical_info['double'] = bonds
            else:
                print("must select only 2 atoms")
                self.selections[selection].selected_atoms[self.click_mode] = []

            self._set_style_monomer_view()
            self._set_clickable('do_nothing')
            self._set_hoverable('do_nothing')
            self._enable_buttons(temp_disable_widgets)
            self._reload_view()
        else:
            print("internal: no action for _button_assign_double_bonds!")
        
    def _function_triple_bonds_button(self, b):
        temp_disable_widgets = ["double_bonds_button", "formal_charge_button", "formal_charge_menu"]
        selection = self._get_selection()
        
        # cycles between "Make Triple bonds" and "Commit Triple Bond"
        if b.description == "Make Triple bonds":
            b.description = "Commit Triple Bond"
            b.button_style = "warning"
            self._set_clickable('select_triple')
            self._set_hoverable('select_triple')
            self.view.center({"serial": self.selections[selection].selected_atoms['select_monomer']})
            self._disable_buttons(temp_disable_widgets)
            self._reload_view()
        elif b.description == "Commit Triple Bond":
            b.description = "Make Triple bonds"
            b.button_style = ""
            
            # load the new triple bond into the view with the chemistry engine + sdf
            # check that only two atoms were selected
            bond = self.selections[selection].selected_atoms['select_triple']
            self.selections[selection].selected_atoms['select_triple'] = []
            if len(bond) == 2:
                new_bond = tuple(bond) # sdf format indexed from 0
                bonds = self.selections[selection].chemical_info['triple'] # does default dict handle this? TODO
                bonds.append(new_bond)
                self.selections[selection].chemical_info['triple'] = bonds
            else:
                print("must select only 2 atoms")
                self.selections[selection].selected_atoms[self.click_mode] = []

            self._set_style_monomer_view()
            self._set_clickable('do_nothing')
            self._set_hoverable('do_nothing')
            self._enable_buttons(temp_disable_widgets)
            self._reload_view()
        else:
            print("internal: no action for _button_assign_triple_bonds!")
        
    def _function_formal_charge_button(self, b):
        print("in function")
        
    def _function_run_button(self, b):
        json_dict = self.chemistry_engine.substructure_generator.get_monomer_info_dict()
        substructure_json = "substructures.json"
        with open(substructure_json, "w") as file:
            json.dump(json_dict, file, indent=4)
        assigned_atoms, _, unassigned_atoms, _, chemical_info, isomorphism_summary = self.chemistry_engine.test_polymer_load(substructure_json)
        if len(unassigned_atoms) == 0:
            self.widgets['valid_check'].value = True
            self.widgets['valid_check'].readout = "Complete"
        else:
            self.widgets['valid_check'].value = False
            self.widgets['valid_check'].readout = "Incomplete"
        self.isomorphisms = isomorphism_summary
        self.highlights = self.highlights | set(unassigned_atoms)
        sdf_block = self.chemistry_engine.get_sdf_block(chemical_info)
        self.view = self._view_from_sdf_block(sdf_block)
        self._add_highlights(self.highlights)
        self._set_clickable('do_nothing')
        self._set_hoverable('do_nothing')
        self.view.setStyle({"model": -1}, {"stick": {}})
        self.widgets['color_code_button'].disabled = False
        self._reload_view()
        
    def _function_print_button(self, b):
        json_dict = self.chemistry_engine.substructure_generator.get_monomer_info_dict()
        print(json.dumps(json_dict, indent=4))

    def _function_color_code_button(self, b):

        if b.button_style == "":
            b.button_style = "info"
            cmap = cm.hsv
            names = set()
            for name, ids, selected in self.isomorphisms:
                names.add(name)
            names = list(names)
            length = len(names)
            inputs = dict()
            for i, name in zip(range(0, length), names):
                inputs[name] = matplotlib.colors.rgb2hex(cmap(i/length))
            self.monomer_colors = inputs
            self.widgets['iso_inspect_dropdown'].options = list(set(self.widgets['iso_inspect_dropdown'].options) | set(names))
            self.widgets['iso_inspect_dropdown'].disabled = False
            self._function_iso_inspect_dropdown({'new': self.widgets["iso_inspect_dropdown"].value}) # update the iso inspect function
        elif b.button_style == "info":
            b.button_style = ""
            self.view.setStyle({"model": -1}, {"stick": {"colorscheme": "default"}})
            self.monomer_colors = dict()
            self._reload_view()

    def _function_remove_highlights_button(self, b):
        self.view.removeAllShapes()
        self.highlights = set()
        self._load_buttons()
        self._reload_view()
           
    def _function_header_tags(self, change):
        if change['new'] == "Edit Monomers":
            sdf_block = self.chemistry_engine.get_sdf_block()
            self.view = self._view_from_sdf_block(sdf_block)
            self._add_highlights(self.highlights)
            self._set_clickable('select_monomer')
            self._set_hoverable('select_monomer')
            self._load_buttons("edit-select_monomer")
            self.view.setStyle({"model": -1}, {"stick": {"colorscheme": "default"}})
        elif change['new'] == "Test Load":
            self._set_clickable('do_nothing')
            self._set_hoverable('do_nothing')
            self._load_buttons("test_load-default")
            
        self._reload_view()
    
    def _function_terminal_group_tags(self, change):
        new_value = change['new']
        if new_value not in self.selections.keys():
            print("error in _function_terminal_group_tags")
            return
        self._set_style_monomer_view()
        self._reload_view()

    def _function_terminal_group_delete_button(self, b):
        current_tag = self._get_selection()
        if current_tag in ["main", "No Terminal Groups Found"]:
            return
        else:
            self.selections.pop(current_tag)
        self._load_buttons()
        self._reload_view()

    def _function_iso_inspect_dropdown(self, change):
        new_value = change['new']
        self.view.setStyle({"model": -1}, {"stick": {"colorscheme": "default"}})
        self.view.removeAllShapes()
        if new_value in self.monomer_colors.keys():
            color = self.monomer_colors[new_value]
            self.widgets['color_box'].value = color
            for name, ids, selected in self.isomorphisms:
                if name == new_value:
                    self.view.setStyle({"serial": ids}, {"stick": {"color": color}})

        elif new_value == 'Fitted Isomorphisms':
            unmapped_atoms = set(range(0, self.chemistry_engine.n_atoms))
            for name, ids, selected in self.isomorphisms:
                if not selected:
                    continue
                unmapped_atoms = unmapped_atoms - set(ids)
                color = self.monomer_colors[name]
                self.view.setStyle({"serial": ids}, {"stick": {"color": color}})
            self._add_highlights(unmapped_atoms)

        elif new_value == 'Unmapped Atoms':
            unmapped_atoms = set(range(0, self.chemistry_engine.n_atoms))
            for name, ids, selected in self.isomorphisms:
                unmapped_atoms = unmapped_atoms - set(ids)
            self._add_highlights(unmapped_atoms)
        else:
            print("internal error")
            return
        self._reload_view()

class Selection:
    def __init__(self):
        self.selected_atoms = defaultdict(list)
        self.chemical_info = defaultdict(list)

if __name__ == "__main__":
    file = "openff_polymer_testing/polymer_examples/rdkit_simple_polymers/dyed_protein.pdb"
    engine = ChemistryEngine(file)
    for bond in engine.full_molecule.GetBonds():
        bond.SetBondType(Chem.BondType.SINGLE)
    viz = PolymerVisualizer3D(engine)
    viz._function_run_button(None)

    # TODO:
    # 1) Better monomer visualization and color-coding
    #   a) use buttons to switch between monomer-specific isomorphisms and "fitting" isomorphisms
    #   b) figure out a way to color code the buttons or have a color map to to the side
    #   c) have a "complete matching" vs "incomplete matching" box somewhere
    # 2) figure out what was up with the error during meeting with Dr. Shirts
    # 3) finalize what to do with caps and how to treat them
    # 4) complete the workflow to end in a simulation with gasteiger charges 
    # 5) [if time], flesh out tool to give library charges