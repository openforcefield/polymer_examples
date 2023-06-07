# Code to populate metadata of openforcefield molecule with partitioned 
# assingment of matches. Must be called after successful call to openff.toolkit.topology.topology.from_pdb()
from openff.toolkit.topology import Topology
import networkx as nx
import itertools
import numpy as np
import json
from collections import defaultdict

def partition(offtop):
    """
    args
    ----
    offtop: Openforcefield toolkit topology with complete metadata following successful call
        to openff.toolkit.topology.topology.from_pdb().
    returns 
    -------
    True if partition was successful
    False otherwise. This does not mean that the resulting assignment is invalid. It only means
        that clean partitioning cannot be achieved and some matches necessarily overlap. Try choosing
        other substructures if partitioning is necessary. 
    """
    
    def _create_choice_graph(isomorphism_info):
        isomorphism_info = list(isomorphism_info.values()) # original ids are of no use
        all_values = set()
        for isomorphism in isomorphism_info:
            iso_ids, adjacencies = isomorphism
            for value in iso_ids:
                all_values.add(value)
        max_val = max(all_values)
        min_val = min(all_values)
        # create a matrix to find atom overlaps
        matrix = np.zeros((len(isomorphism_info), max_val+1))
        adjacency_list = []
        for idx, isomorphism in enumerate(isomorphism_info):
            iso_ids, adjacencies = isomorphism
            cap_ids = []
            adjacencies_tupled = [] # adjacencies in tuple form for easier searching
            for bond, bond_type in adjacencies.items():
                edge_atom, cap_atom = bond
                cap_ids.append(cap_atom)
                adjacencies_tupled.append(tuple([edge_atom, cap_atom, bond_type]))
            adjacency_list.append(adjacencies_tupled)

            matrix[idx, list(set(iso_ids)-set(cap_ids))] = 1

        # matrix_data and adjacency_list are parallel to eachother
        # each item in adjacency_list corresponds to the isomorphism found in the 
        # corresponding row of matrix_data
        n_rows, n_cols = matrix.shape
        G = nx.Graph()
        all_bonds = []
        for iso_id, isomorphism in enumerate(adjacency_list):
            n_atoms = len(isomorphism_info[iso_id][0])
            if isomorphism == []:
                continue
            neighbors = []
            for bond in isomorphism:
                edge_a, cap_a, bond_type = bond
                all_bonds.append((edge_a, cap_a, bond_type, iso_id))

            row, = np.nonzero(np.squeeze(np.asarray(matrix[iso_id, :])))
            row = list(row)
            # find overlaps using matrix 
            overlaps, = np.nonzero(np.asarray(matrix[:, row].sum(axis=1)).reshape(-1))
            overlaps = list(overlaps)
            overlaps.remove(iso_id)
            G.add_node(
                iso_id,
                selected = False,
                overlaps = overlaps,
                bonds = isomorphism,
                n_atoms = n_atoms
            )
        for bond1, bond2 in itertools.combinations(all_bonds, 2):
            edge_a1, cap_a1, bond_type1, iso_id1 = bond1
            edge_a2, cap_a2, bond_type2, iso_id2 = bond2
            if (edge_a2, cap_a2, bond_type2) == (cap_a1, edge_a1, bond_type1):
                G.add_edge(
                    iso_id1,
                    iso_id2,
                    order = bond_type1
                )

        return G

    def _traverse(starting_queue, graph):
        found_nodes = []
        found_edges = []
        queue = starting_queue
        while len(queue) > 0:
            if queue == []:
                return found_nodes
            v = queue.pop(0)
            found_nodes.append(v)
            for bond in graph.nodes[v]['bonds']:
                edge_a, cap_a, bond_type = bond
                if {edge_a, cap_a} in found_edges:
                    continue
                selected_neighbor = -1
                for neighbor in graph.neighbors(v):
                    # make sure this is the neighbor with the correct bond info
                    neighbor_node = graph.nodes[neighbor]
                    if (cap_a, edge_a, bond_type) not in neighbor_node['bonds']:
                        continue
                    if neighbor_node['selected'] == True:
                        selected_neighbor = -1
                        break
                    selected_neighbor = neighbor
                if selected_neighbor >= 0:
                    queue.append(selected_neighbor)
                    found_edges.append({edge_a, cap_a})

        return found_nodes

    for offmol in offtop.molecules:
        # first, populate isomorphism_info, which has the following format:
        #  [([target_ids], [{(map_start, map_end): bond_type}])]
        # isomorphism_info = dict()
        target_ids = defaultdict(list)
        bond_info = defaultdict(dict)
        match_names = defaultdict(str)
        for atom in offmol.atoms:
            iso_info = {int(k): entry for k,entry in json.loads(atom.metadata["match_info"]).items()}
            for query_num, query_name_id in iso_info.items():
                match_name = query_name_id[0] # id of the query atom
                query_id = query_name_id[1] # total number of query 
                match_names[query_num] = match_name
                if atom.molecule_atom_index not in target_ids[query_num]:
                    target_ids[query_num].append(atom.molecule_atom_index)
                for b in atom.bonds:
                    # try to find inter-monomer bonds
                    if b.atom1_index == atom.molecule_atom_index:
                        begin_atom = b.atom1
                        end_atom = b.atom2
                    else:
                        begin_atom = b.atom2
                        end_atom = b.atom1
                    if query_num in [int(k) for k in json.loads(end_atom.metadata["match_info"]).keys()]:
                        neighbor = False
                    else:
                        neighbor = True
        
                    if neighbor:
                        bond_entry = tuple([begin_atom.molecule_atom_index, end_atom.molecule_atom_index])
                        if bond_entry not in bond_info[query_num]:
                            bond_info[query_num][bond_entry] = b.bond_order

        assert len(target_ids) == len(bond_info)
        isomorphism_info_dict = {idx: (target_ids[idx], bond_info[idx]) for idx in target_ids}
        idx_to_match_id = {idx: match_id for idx, match_id in enumerate(target_ids)}
        choice_G = _create_choice_graph(isomorphism_info_dict)
        connected_component_counts = []
        biggest_chain = None
        biggest_chain_length = 0
        for chain in nx.connected_components(choice_G):
            subgraph = choice_G.subgraph(chain)
            not_searched_nodes = list(subgraph.nodes)
            while len(not_searched_nodes) != 0:
                unique_group = _traverse([not_searched_nodes[0]], subgraph)
                new_tally = 0
                old_tally = 0
                overlapping_nodes = []
                for node in unique_group:
                    new_tally += subgraph.nodes[node]['n_atoms']
                    for overlapping_node in subgraph.nodes[node]['overlaps']:
                        if overlapping_node in subgraph.nodes and subgraph.nodes[overlapping_node]['selected'] == True:
                            overlapping_nodes.append(overlapping_node)
                            old_tally += subgraph.nodes[overlapping_node]['n_atoms']
                if new_tally > old_tally:
                    # exchange new for old choices 
                    for node in unique_group:
                        subgraph.nodes[node]['selected'] = True
                    for node in overlapping_nodes:
                        subgraph.nodes[node]['selected'] = False

                [not_searched_nodes.remove(i) for i in unique_group if i in not_searched_nodes]

            tally = 0
            for node in subgraph.nodes:
                if subgraph.nodes[node]['selected']:
                    tally += subgraph.nodes[node]['n_atoms']
            if tally > biggest_chain_length:
                biggest_chain = subgraph
                biggest_chain_length = tally
            # finally pick from unique_mapping_groups the largest list
        largest_mapping_group = []
        for node in biggest_chain.nodes:
            if biggest_chain.nodes[node]['selected']:
                largest_mapping_group.append(node)

        all_ids = []
        for iso_id in largest_mapping_group:
            all_ids += target_ids[idx_to_match_id[iso_id]]

        if set(all_ids) == set([a.molecule_atom_index for a in offmol.atoms]):
            # assign metadata
            for iso_id in largest_mapping_group:
                res_id = idx_to_match_id[iso_id]
                molecule_ids = target_ids[res_id]
                match_name = match_names[res_id]
                for m_id in molecule_ids:
                    atom = offmol.atom(m_id)
                    iso_info = {int(k): entry for k,entry in json.loads(atom.metadata["match_info"]).items()}
                    picked_match = iso_info[res_id]
                    # new metadata
                    atom.metadata["residue_name"] = picked_match[0]
                    atom.metadata["substructure_query_id"] = picked_match[1]
                    atom.metadata["residue_number"] = res_id
        else:
            return False # if one mol fails, return False
    return True # if all are successful, return True