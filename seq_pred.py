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
from collections import Counter
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
             
inv_sets = list(set() for i in neighbors)            
inverse = dict(zip(list(neighbors.keys()),inv_sets)) # inverse tells us for each residue, which TERMs include it
        
for i in neighbors:
    neighbors[i] = set(neighbors[i]) 

overlap = copy.deepcopy(neighbors) # overlap tells us, for each TERM, which TERMs have shared residues with it

# For each residue, find which TERMs include it, all these TERMs should be connected to each other           
for i in neighbors:
    for j in neighbors[i]:
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
    def __init__(self,seq):
        self.seq = seq
        
    def select(self,node,number):
        return self.seq[node][number]
    
    def nb_frag(self, node):
        return len(self.seq[node])
    
for i in neighbors: 
    G.add_node(i, matches = node_attributes(match_sequence).seq[i])
     

# add edges
edges = []
for i in overlap:
    for j in overlap[i]:
        edges.append((i,j))

G.add_edges_from(edges)       
G.remove_edges_from(G.selfloop_edges())   

# add edge attributes
for i in G.edges():
    G.add_edge(i[0],i[1],sameAA = find_overlap_position(neighbors_copy[i[0]],neighbors_copy[i[1]]))


for i in range(10):
    pop = geneticAlgorithm(100, 50, 2, 0.03, node_attributes(match_sequence), 20, G)
    np.savetxt('result' + str(i) + '.txt',pop[0:3,:], fmt = '%i')


#geneticAlgorithmPlot(100, 50, 2, 0.03, node_attributes(match_sequence), 500, G)

possible_res = list([] for i in neighbors)            
candidate = dict(zip(list(neighbors.keys()),possible_res)) # inverse tells us for each residue, which TERMs include it

predict = dict(zip(keys,list(pop[0]))) # predict represents the choice of fragment for each TERM
for i in inverse:
    for j in inverse[i]:
        place = neighbors[j].index(i)
        nb_fragment = predict[i]
        if nb_fragment >= len(node_attributes(match_sequence).seq[j]):
            print(j + ' is a gap')
            continue
        else:
            candidate[i].append(node_attributes(match_sequence).select(j,nb_fragment)[place])


possible_seq = ''
for i in keys:
    if candidate[i] != []:
        possible_seq += Counter(candidate[i]).most_common(1)[0][0]
    else:
        possible_seq += '-'
        continue
    
