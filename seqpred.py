#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 11 11:40:57 2019

@author: pengdandan
"""
import copy, glob, re
import numpy as np
import networkx as nx

from prody import *
from Bio.PDB.Polypeptide import three_to_one


## sort string based on the embedded number

def embedded_numbers(s):
    re_digits = re.compile(r'(\d+)')
    pieces = re_digits.split(s)                 
    pieces[1::2] = map(int, pieces[1::2])       
    return pieces
 
def sort_string(lst):
    return sorted(lst, key=embedded_numbers)    

## find corresponding amino acids for u, v on edge (u,v)
def find_overlap_position(term1,term2):
    position = []
    for i in set(term1) & set (term2):
        position.append((term1.index(i),term2.index(i)))
    return position


# read TERMs and find structure overlap between TERMs

def get_neighbors(protein_id):
    #path = '/cluster/home/pengd/project/train_data/'+ protein_id + '/fragments' # path of folder where stores the TERMs' PDB file 
    path = '/Users/pengdandan/Desktop/lab_rotation/LabRotation2/data/' + protein_id
    protein = parsePDB(protein_id + '.pdb').select('protein').copy()
    residues = []
    neighbors = {} #store the amino acids that consitute each TERM
    for res in protein.iterResidues():
        residues.append(res)
        cid, resnum = res.getChid(),res.getResnum()
        neighbors[cid + str(resnum)] = []
        term_frag = path + '/' + protein_id + '_' + cid + str(resnum) + '.pdb'
        f_struct = parsePDB(term_frag)
        for fres in f_struct.iterResidues():
            fcid, fresnum = fres.getChid(),fres.getResnum()
            neighbors[cid + str(resnum)].append(fcid + str(fresnum))
            
    return neighbors

            

def get_inv_and_overlap(keys, neighbors):
    empty_sets = list(set() for i in neighbors)            
    inverse = dict(zip(list(keys),empty_sets)) # inverse tells us for each residue, which TERMs include it

    neigh_sets = list(set(neighbors[i]) for i in keys)     
    neighbors_sets = dict(zip(list(keys), neigh_sets))

    for i in neighbors:    
        neighbors[i] = sort_string(neighbors[i]) 

    overlap = copy.deepcopy(neighbors_sets) # overlap tells us, for each TERM, which TERMs have shared residues with it

# For each residue, find which TERMs include it, all these TERMs should be connected to each other           
    for i in neighbors_sets:
        for j in neighbors_sets[i]:
            inverse[j].add(i) 
      
    for i in inverse:
        for j in inverse[i]:
            for k in inverse[i] - {j}:
                overlap[j].add(k)    
    
    return inverse, overlap
                  
# read match sequences
def get_match(protein_id):
    #dsScore_uniq = glob.glob(r'/cluster/home/pengd/project/train_data/' + protein_id + '/designscore/uniq*.seq')     
    dsScore_uniq = glob.glob(r'/Users/pengdandan/Desktop/lab_rotation/LabRotation2/code/2QMT_designscore/uniq*.seq')
    dsScore_uniq.sort() 

    match = {}

    for file in dsScore_uniq:
        name = (file.split('_')[-1]).split('.')[0]
        match[name] = np.loadtxt(file, dtype = np.str)
        
        if len(match[name].shape) > 1:
            match[name] = match[name][:,1:]
        else:
            match[name] = match[name][1:]

    # convert unnatural amino acids to common ones
    # convert three-letter amino acids to one-letter ones
    toCommonAA = {'ASX':'ASN','CSO':'CYS','GLX':'GLU','HIP':'HIS','HSC':'HIS',
              'HSD':'HIS','HSE':'HIS','HSP':'HIS','MSE':'MET','SEC':'CYS',
              'SEP':'SER','TPO':'THR','PTR':'TYR','XLE':'LEU'}

    for i in match:
        for j in toCommonAA:
            match[i] = np.where(match[i] == j,toCommonAA[j], match[i])
        match[i] = np.vectorize(three_to_one)(match[i])   
    
    return match


class node_attributes:
    def __init__(self, seq):
        self.seq = seq
        
    def select(self, node, number):
        if len(self.seq[node].shape) > 1:
            return self.seq[node][number]
        else:
            return self.seq[node]
            
    def count(self, node):
        if len(self.seq[node].shape) > 1:
            return len(self.seq[node])
        else:
            return 1
  
# build graph for connectivity between TERMs
def build_graph(overlap, frag, neighbors):
    G = nx.Graph()
    # add nodes with attributes
    for i in neighbors: 
        G.add_node(i, match = frag.seq[i])
    
    # add edges
    edges = []
    for i in overlap:
        for j in overlap[i]:
            edges.append((i,j))

    G.add_edges_from(edges)       
    G.remove_edges_from(G.selfloop_edges())  
    
    # add edge attributes
    for i in G.edges():
        G.add_edge(i[0],i[1],sameAA = find_overlap_position(neighbors[i[0]],neighbors[i[1]]))

    return G

def get_information(protein_id):
    neighbors = get_neighbors(protein_id)
    keys = sort_string(neighbors.keys()) 
    inverse, overlap = get_inv_and_overlap(keys, neighbors)
    match = get_match(protein_id)

    frag = node_attributes(match)   

    count = []
    for i in keys:
        count.append(frag.count(i))
        frags_count = dict(zip(keys, count))

    G = build_graph(overlap, frag, neighbors)
    
    return neighbors, keys, inverse, overlap, match, frag, frags_count, G
