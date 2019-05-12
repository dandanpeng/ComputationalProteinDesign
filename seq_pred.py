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
from random import sample, choice
from Bio import Align
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
    # protein_pdb_file: path of the protein PDB file
    # path_of_TERMs_file: path of the folder which include TERMs file
    # protein_id: ID of the protein
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

for i in node_attributes:
    node_attributes[i] = np.where(node_attributes[i] == 'ASX',
                                  'ASN', node_attributes[i])
    node_attributes[i] = np.where(node_attributes[i] == 'CSO',
                                  'CYS', node_attributes[i])
    node_attributes[i] = np.where(node_attributes[i] == 'GLX',
                                  'GLU', node_attributes[i])
    node_attributes[i] = np.where(node_attributes[i] == 'HIP',
                                  'HIS', node_attributes[i])
    node_attributes[i] = np.where(node_attributes[i] == 'HSC',
                                  'HIS', node_attributes[i])
    node_attributes[i] = np.where(node_attributes[i] == 'HSD',
                                  'HIS', node_attributes[i])
    node_attributes[i] = np.where(node_attributes[i] == 'HSE',
                                  'HIS', node_attributes[i])
    node_attributes[i] = np.where(node_attributes[i] == 'HSP',
                                  'HIS', node_attributes[i])
    node_attributes[i] = np.where(node_attributes[i] == 'MSE',
                                  'MET', node_attributes[i])
    node_attributes[i] = np.where(node_attributes[i] == 'SEC',
                                  'CYS', node_attributes[i])
    node_attributes[i] = np.where(node_attributes[i] == 'SEP',
                                  'SER', node_attributes[i])
    node_attributes[i] = np.where(node_attributes[i] == 'TPO',
                                  'THR', node_attributes[i])
    node_attributes[i] = np.where(node_attributes[i] == 'PTR',
                                  'TYR', node_attributes[i])
    node_attributes[i] = np.where(node_attributes[i] == 'XLE',
                                  'LEU', node_attributes[i])
for i in node_attributes:
    for x in range(node_attributes[i].shape[0]):
        for y in range(node_attributes[i].shape[1]):
            node_attributes[i][x,y] = three_to_one(node_attributes[i][x,y])

# add edge with attributes
edges = []
for i in overlap:
    for j in overlap[i]:
        edges.append((i,j))

G.add_edges_from(edges)       
G.remove_edges_from(G.selfloop_edges())   

for i in G.edges:
    G.add_edge(i[0],i[1],sameAA = find_overlap_position(neighbors_copy[i[0]],neighbors_copy[i[1]]))
 
 


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
        
        
def crossover(parent1, parent2,num_points):
    points = sample([i for i in range(len(parent1))], num_points)
    points = list(set(points) | {0,len(parent1)-1})
    points.sort()
    child = []
    i = 0  
    names = locals()
    while i < len(points)-1: 
        names['child' + str(i)] = choice([parent1[points[i]:points[i+1]],parent2[points[i]:points[i+1]]])
        child = child + names.get('child'+str(i))
        i+=1  
    return child
     
        
def mutate(gene,mutationRate):
    mutationNumber = int(mutationRate * len(gene))
    mutationPosition = sample(gene,mutationNumber)
    for i in mutationPosition:
        gene[i] = choice(len(node_attributes['A'+str(i)]))
    return gene


def get_fitness(individual):
    aligner = Align.PairwiseAligner()
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    aligner.substitution_matrix = blosum62
    for i in range(individual):
        
        
        
    return aligner.score(seq1,seq2)

    

def checkOverlap(s1,s2):
    m = min(len(s1),len(s2))
    for i in range(m,0,-1):
        if s1[-i:] == s2[:i]:
            #print('TRUE')
            #return s1+s2[i:]
            return True
