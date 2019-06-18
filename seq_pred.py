#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 11 11:40:57 2019

@author: pengdandan
"""
import copy, glob
import pandas as pd, numpy as np
import networkx as nx

from prody import *
from function import *
from Bio.PDB.Polypeptide import three_to_one

# read TERMs and find structure overlap between TERMs
path = '/cluster/home/pengd/project/test/2QMT/fragments' # path of folder where stores the TERMs' PDB file 
    
protein_id = '2QMT'
protein = parsePDB('2QMT.pdb').select('protein').copy()
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

keys = sort_string(neighbors.keys())             

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

# neighbors_copy stores the amino acids which constitute each TERM, 
# but the amino acid which is same with the one who defined the TERM is removed
neighbors_copy = copy.deepcopy(neighbors)

for i in neighbors_copy:
    neighbors_copy[i] = sort_string(neighbors_copy[i])
    neighbors_copy[i].remove(i)
    
                       
# read match sequences
dsScore_uniq = glob.glob(r'/cluster/home/pengd/project/test/2QMT/designscore/uniq*.seq')     
dsScore_uniq = glob.glob(r'/Users/pengdandan/Desktop/lab_rotation/LabRotation2/test/2QMT /designscore/uniq*.seq')
dsScore_uniq.sort() 

match_sequence = {}

for file in dsScore_uniq:
    name = (file.split('_')[-1]).split('.')[0]
    df = pd.read_csv(file,header = None, sep = '\t')
    seq = pd.DataFrame(df.iloc[:,0].str.split().tolist()).iloc[:,1:]
    match_sequence[name] = seq.values

# convert unnatural amino acids to common ones
toCommonAA = {'ASX':'ASN','CSO':'CYS','GLX':'GLU','HIP':'HIS','HSC':'HIS',
              'HSD':'HIS','HSE':'HIS','HSP':'HIS','MSE':'MET','SEC':'CYS',
              'SEP':'SER','TPO':'THR','PTR':'TYR','XLE':'LEU'}

for i in match_sequence:
    for j in toCommonAA:
        match_sequence[i] = np.where(match_sequence[i] == j,toCommonAA[j],match_sequence[i])

# convert three-letter amino acids to one-letter ones
for i in match_sequence:
    for x in range(match_sequence[i].shape[0]):
        for y in range(match_sequence[i].shape[1]):
            match_sequence[i][x,y] = three_to_one(match_sequence[i][x,y])

# build graph for connectivity between TERMs
G = nx.Graph()

# add node with attributes
class node_attributes:
    def __init__(self, seq):
        self.seq = seq
        
    def select(self, node, number):
        return self.seq[node][number]
    
    def count(self, node):
        return len(self.seq[node])

    
for i in neighbors: 
    G.add_node(i, matches = fragments.seq[i])
     

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


#sequence = np.zeros(shape = (30, 56), dtype = int)
#for i in range(10):
    #pop = geneticAlgorithm(100, 50, 2, 0.03, node_attributes(match_sequence), 50, G)
    #sequence[0 + 3*i :3 + 3*i] = pop[0:3]

#f = open("/Users/pengdandan/Desktop/lab_rotation/LabRotation2/test/sequence2.txt",'w')

#for i in range(len(sequence)):
    #f.write('>seq' + str(i) + '\n')
    #f.write(restore_seq(keys, sequence[i], neighbors,inverse, node_attributes(match_sequence)) + '\n')
    
#f.close()

plot(100, 50, 3, 0.03, node_attributes(match_sequence), 50, G)
plot(100, 50, 3, 0.03, node_attributes(match_sequence), 500, G)

