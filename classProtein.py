#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 11 11:40:57 2019

@author: pengdandan
"""
import re
import numpy as np
import networkx as nx

from prody import *
from Bio.PDB.Polypeptide import three_to_one

path = path = '/Users/pengdandan/Desktop/lab_rotation/LabRotation2/data/2QMT'
match_path = '/Users/pengdandan/Desktop/lab_rotation/LabRotation2/code/2QMT_designscore/uniq_t1k_'

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

toCommonAA = {'ASX':'ASN','CSO':'CYS','GLX':'GLU','HIP':'HIS','HSC':'HIS',
              'HSD':'HIS','HSE':'HIS','HSP':'HIS','MSE':'MET','SEC':'CYS',
              'SEP':'SER','TPO':'THR','PTR':'TYR','XLE':'LEU'}

class Protein:
    
    def __init__(self, protein_id):
        self.protein_id = protein_id    
        self.neighbors = None
        self.terms = None
        self.inverse = None
        self.overlap = None
        self.match = None
        self.graph = None
        
        self.get_neighbors()
        self.get_term_name()
        self.get_inv_and_overlap()
        self.get_match()
        self.build_graph()
    # read TERMs & store in dict 'neighbors'
    def get_neighbors(self):
        #path = '/cluster/home/pengd/project/train_data/'+ self.protein_id + '/fragments' # path of folder where stores the TERMs' PDB file 
        protein = parsePDB(self.protein_id + '.pdb').select('protein').copy()
        neighbors = {} #store the amino acids that consitute each TERM
        for res in protein.iterResidues():
            cid, resnum = res.getChid(),res.getResnum()
            neighbors[cid + str(resnum)] = []
            term_frag = path + '/' + self.protein_id + '_' + cid + str(resnum) + '.pdb'
            f_struct = parsePDB(term_frag)
            for fres in f_struct.iterResidues():
                fcid, fresnum = fres.getChid(),fres.getResnum()
                neighbors[cid + str(resnum)].append(fcid + str(fresnum))
        self.neighbors = neighbors
        
    def get_term_name(self):
        terms = []
        for i in self.neighbors:
            terms.append(i)
        self.terms = terms
               
    def get_inv_and_overlap(self):
        empty_sets = list(set() for i in self.terms)            
        inverse = dict(zip(list(self.terms), empty_sets)) # inverse tells us for each residue, which TERMs include it
        #neigh_sets = list(set(self.neighbors[i]) for i in self.terms)     
        #neighbors_sets = dict(zip(list(self.terms), neigh_sets))

        #for i in self.neighbors:    
            #self.neighbors[i] = sort_string(self.neighbors[i]) 
        neighbors_set = self.neighbors.copy() # overlap tells us, for each TERM, which TERMs have shared residues with it
        for i in self.terms:
            neighbors_set[i] = set(self.neighbors[i])
        overlap = neighbors_set.copy()
        
        # For each residue, find which TERMs include it, all these TERMs should be connected to each other           
        for i in neighbors_set:
            for j in neighbors_set[i]:
                inverse[j].add(i)       
        for i in inverse:
            for j in inverse[i]:
                for k in inverse[i] - {j}:
                    overlap[j].add(k)    
    
        self.inverse = inverse
        self.overlap = overlap

    # read match sequences
    def get_match(self): 
        #match_path = '/cluster/home/pengd/project/train_data/' + self.protein_id + '/designscore/uniq_t1k_' + self.protein_id     
        match = {}
        for i in range(len(self.terms)):
            term = self.terms[i]
            file = match_path + self.protein_id + '_' + term +'.seq'
            match[term] = np.loadtxt(file, dtype = np.str)        
            if len(match[term].shape) > 1:
                match[term] = match[term][:,1:]
            else:
                match[term] = match[term][1:]
        # convert unnatural amino acids to common ones
        # convert three-letter amino acids to one-letter ones        
        for i in self.terms:
            for j in toCommonAA:
                match[i] = np.where(match[i] == j,toCommonAA[j], match[i])
            match[i] = np.vectorize(three_to_one)(match[i])   
    
        self.match = match

    def build_graph(self):
        G = nx.Graph()
        '''
        # add nodes with attributes
        for i in self.terms: 
            G.add_node(i, match = self.match[i])
        '''    
        # add edges
        edges = []
        for i in self.overlap:
            for j in self.overlap[i]:
                edges.append((i,j))

        G.add_edges_from(edges)       
        G.remove_edges_from(G.selfloop_edges())  
    
        # add edge attributes
        for i in G.edges():
            G.add_edge(i[0],i[1],sameAA = find_overlap_position(self.neighbors[i[0]],self.neighbors[i[1]]))

        self.graph = G
    

class Matches:
    
    def __init__(self, seq):
        self.seq = seq
        self.fargs_count = None
        
        self.count_frags()

    def count_frags(self):
        count = {}
        for i in self.seq:
            if len(self.seq[i].shape) > 1:
                count[i] = len(self.seq[i])
            else:
                count[i] = 1
        self.frags_count = count
        
    def select_frag(self, term, number):
        if number < self.frags_count[term]:
            if len(self.seq[term].shape) > 1:
                return self.seq[term][number]
            else:
                return self.seq[term]
        else:
            return None
        

