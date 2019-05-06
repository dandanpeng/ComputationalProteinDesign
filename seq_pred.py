#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 11 11:40:57 2019

@author: pengdandan
"""
import copy
import re
import pandas as pd
import glob
import numpy as np
import networkx as nx

from prody import *
from Bio import pairwise2
from Bio.SubsMat.MatrixInfo import blosum62
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder,three_to_one

##### read original sequence
structure = PDBParser().get_structure('2QMT','2QMT.pdb')
ppb = PPBuilder()
for pp in ppb.build_peptides(structure):
    ref = str(pp.get_sequence())
            
# read TERMs
path = r'/Users/pengdandan/Desktop/lab_rotation/LabRotation2/test/2QMT/'
all_files =glob.glob("*.pdb")
all_files.sort()

terms = []
for filename in all_files:
    terms.append((filename.split('.')[0]).split('_')[1])
        
# find structure overlap between TERMs
def checkTermOverlap(protein_pdb_file,path_of_TERMs_file,protein_id):
    protein = parsePDB(protein_pdb_file).select('protein').copy()
    residues = []
    neighbors = {}
    for res in protein.iterResidues():
        residues.append(res)
        cid,resnum = res.getChid(),res.getResnum()
        neighbors[cid+','+str(resnum)] = []
        term_frag = path_of_TERMs_file + '/' + protein_id + '_' + cid + str(resnum) + '.pdb'
        f_struct = parsePDB(term_frag)
        for fres in f_struct.iterResidues():
            fcid, fresnum = fres.getChid(),fres.getResnum()
            neighbors[cid + ',' + str(resnum)].append(fcid + ',' + str(fresnum))
            
    keys = list(neighbors.keys())    
    #for i in range(len(keys)-1):
        #already_connect = neighbors[keys[i]]
        #uncompared = set(keys[i+1:]) - set(already_connect) - set(keys[i])
        #for j in list(uncompared):
            #if set(already_connect).intersection(set(neighbors[j])) != set():
                #overlap[keys[i]].append(j)
                #overlap[j].append(keys[i])        
    inv_sets = list(set() for i in keys)            
    inverse = dict(zip(keys,inv_sets))
        
    for i in keys:
        neighbors[i] = set(neighbors[i]) 
    overlap = copy.deepcopy(neighbors)
    # For each residue, find which TERMs include it, all these TERMs should be connected to each other           
    for i in keys:
        for j in neighbors[i]:
            inverse[j].add(i) 
        
    for i in inverse:
        for j in inverse[i]:
            for k in inverse[i] - {j}:
                overlap[j].add(k)

## sort the residue based on numbers in string
re_digits = re.compile(r'(\d+)')
 
def embedded_numbers(s):
    pieces = re_digits.split(s)                 
    pieces[1::2] = map(int, pieces[1::2])       
    return pieces
 
def sort_string(lst):
    return sorted(lst, key=embedded_numbers)    

neighbors_copy = copy.deepcopy(neighbors)
for i in neighbors_copy:
    neighbors_copy[i] = sort_string(neighbors_copy[i])
    neighbors_copy[i].remove(i)
    
               
def find_overlap_position(term1,term2):
    position = []
    for i in set(term1) & set (term2):
        position.append((term1.index(i),term2.index(i)))
    return position
                    
                
        
# read match fragments  
dsScore_uniq = glob.glob('./designscore/uniq*.seq')     
dsScore_uniq.sort() 

# build graph for connectivity between TERMs
G = nx.Graph()
node_attributes = {}

for file in dsScore_uniq:
    name = (file.split('_')[-1]).split('.')[0]
    df = pd.read_table(file,header = None)
    seq = pd.DataFrame(df.iloc[:,0].str.split().tolist()).iloc[:,1:]
    G.add_nodes_from([name],match = seq.values)
    node_attributes[name] = seq.values
    
# add edge with attributes
for i in overlap:
    for j in overlap[i]:
        G.add_edge(i,j,sameAA = find_overlap_position(neighbors_copy[i],neighbors_copy[j]))
 
G.remove_edges_from(G.selfloop_edges())    


# randomly generate sequences as parent
def InitParents(node_attributes,num_parents,len_gene):
    pop = np.zeros(shape = (num_parents,len_gene))
    
    for i in range(num_parents):
        gene = []
        for j in keys:
            gene.append(round(np.random.random(1)[0] * len(node_attributes[j])))
        pop[i] = gene
    return pop

def select(pop, fitness):
    pass
        
        
def crossover(parent1, parent2):
       
        pass
    
def mutate(child):
    pass


def get_fitness(ref,sample):
    return pairwise2.align.localds(ref,sample,blosum62, -10, -1,score_only = True)
    

def checkOverlap(s1,s2):
    m = min(len(s1),len(s2))
    for i in range(m,0,-1):
        if s1[-i:] == s2[:i]:
            #print('TRUE')
            #return s1+s2[i:]
            return True
