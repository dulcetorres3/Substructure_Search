import networkx as nx
from provided import *
import numpy as np
import hashlib
import random
import argparse
import os
from abc import ABC, abstractmethod
import matplotlib.pyplot as plt
from IPython.display import display

class MoleculeGraph:

    color_dict = {
        'C' : 'black',
        'H' : 'white',
        'N' : 'blue',
        'O' : 'red',
        'S' : 'yellow'
    }
    def __init__(self, nodes, bonds, included_hydrogen: Optional[bool] = False):
        if not isinstance(nodes, list):
            raise TypeError("nodes must be in a list")
        
        if not isinstance(bonds, list):
            raise TypeError("bonds must be in a list")
        
        self.hydrogens = included_hydrogen
        self.nodes = nodes
        self.bonds = bonds
        self.graph = nx.Graph()
        self.graph.add_nodes_from(self.nodes)
        self.graph.add_weighted_edges_from(self.bonds)
        self.hetero_atoms = ['S', 'N', 'O']

    @property
    def visualization(self):
        """
        Visualizes molecule instance according to CPK coloring.

        Parameters:
        None

        Returns:
        Visual (matplotlib plot)
        
        """
        colors = [self.color_dict[node[0]] for node in self.nodes]
        visual = nx.draw(self.graph, with_labels=True, node_color=colors)

        display(visual)

    

    @classmethod
    def from_sdf(cls, filename : str, include_hydrogen: Optional[bool] = False):

        full_path = os.path.join('test_compounds', filename)
        try:
            open(os.path.join('test_compounds', filename))
            nodes, bonds = parse_sdf(full_path, include_hydrogen=include_hydrogen)
            nodes = list(nodes.keys())
            return cls(nodes, bonds, included_hydrogen=include_hydrogen)
        except FileNotFoundError:
            print("File is not in the test_compounds sub directory.")

class MoleculeScreen(MoleculeGraph):

    
    def fingerprint(self):

        size = 1000
        cutoff = 6 
        fingerprint = [0] * size

        for node_1 in self.nodes:

            for node_2 in self.nodes:
      

                for path in sorted(nx.all_simple_edge_paths(self.graph, node_1, node_2, cutoff=cutoff)):
                    path_string = " "

                    for node_a, node_b in path:
                        path_string+= node_a[0] + node_b[0] + str(self.graph[node_a][node_b]['weight'])

                    hash_value = hash(path_string)
                    random.seed(hash_value)
                    random_bits = [random.randint(0, 1000) for _ in range(3)]

                    for bit in random_bits:
                            index = bit % size 
                            #index = bit 
                            fingerprint[index] = 1
        return fingerprint


    def sp2_hybridized(self, node:str):
        
        outgoing_edges = self.graph.edges(node, data=True)
        double_bond = False


        for u, v, data in outgoing_edges:
        
            weight = data['weight']
            if weight==2:
                double_bond = True
        
        return double_bond
    
    def hetero_sp2_hybridized(self, node:str):

        # check outgoing edge weights
        outgoing_edges = self.graph.edges(node, data=True)
        double_bond = False
        neighboring_nodes = []

        for u, v, data in outgoing_edges:
            neighboring_nodes.append(v)
            
            weight = data['weight']
            if weight==2:
                double_bond = True

        # check potential delocalization
        if double_bond == False:
            for node in neighboring_nodes:
                double_bond = self.sp2_hybridized()
        
        return double_bond


    def huckel_electrons(self, electrons):
        print(f'electrons: {electrons}')
        test = (electrons - 2) % 4
    
        if test == 0:
            return True
        else:
            return False
        
    
    def aromatic_ring_detection(self):
        
        aromatic_rings = 0
        rings = nx.cycle_basis(self.graph)
        
        pi_electrons = 0
        for ring in rings:
            
            for node in ring:
            
                if node in self.hetero_atoms:
                    pi_electrons += self.hetero_sp2_hybridized(node)
                else:
                    pi_electrons += self.sp2_hybridized(node)
                
            aromatic_rings+= self.huckel_electrons(pi_electrons)


        if aromatic_rings == len(rings):
            return True
        else:
            return False
    
    def screen(self, substructure):

        substructure_fingerprint = np.array(substructure.fingerprint())
        compound_fingerprint = np.array(self.fingerprint())
        
        screen_test = compound_fingerprint - substructure_fingerprint
        
        return np.all(screen_test >= 0)
    
  
        
        

if __name__ == "__main__": 
    
    parser = argparse.ArgumentParser(description= "This script can be used to perform a substructure screen for compounds in the test_compounds subdirectory.")
    parser.add_argument("-sdf_file", help="The name of the file correspnidng to a compound and a substructure in the test_compounds sub directory", type=str, nargs=2)
    parser.add_argument("-nodes_and_bonds", help="A list corresponding to the nodes and bonds of a compound and substructure. Put the list within quotations to avoid error from the commnad line", nargs=4, type=list)
    parser.add_argument("--hydrogens", help="A boolean specifying whether or not hydrigens are included in the compound and substructure", default=False, type=bool)

    args = parser.parse_args()
    if args.sdf_file:
        compound_sdf = args.sdf_file[0]
        substructure_sdf = args.sdf_file[1]
        if args.hydrogens:
            compound = MoleculeScreen.from_sdf(compound_sdf, True)
            substructure = MoleculeScreen.from_sdf(substructure_sdf, True)
        else:
            compound = MoleculeScreen.from_sdf(compound_sdf, True)
            substructure = MoleculeScreen.from_sdf(substructure_sdf)

    elif args.nodes_and_bonds:
        compound_nodes = args.nodes_and_bonds[0]
        compound_bonds = args.nodes_and_bonds[1]
        substructure_nodes = args.nodes_and_bonds[2]
        substructure_bonds = args.nodes_and_bonds[3]
        if args.hydrogens:
            compound = MoleculeScreen(compound_nodes, compound_bonds, included_hydrogen=True)
            substructure = MoleculeScreen(substructure_nodes, substructure_bonds, included_hydrogen=True)
        else:
            compound = MoleculeScreen(compound_nodes, compound_bonds)
            substructure = MoleculeScreen(substructure_nodes, substructure_bonds, substructure_bonds)


    screen_results = compound.screen(substructure)
    if screen_results ==True:
        print("The substructure is present in the compound.")
    else:
        print("The substructure is not present in the compound.")