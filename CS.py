#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  9 16:24:37 2019

@author: pengdandan
"""
import math
import numpy as np

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

from Bio import Align
from Bio.PDB.PDBParser import PDBParser
from Bio.SubsMat.MatrixInfo import blosum62
from Bio.PDB.Polypeptide import PPBuilder

from classProtein import *
from collections import Counter
from random import random


def cuckoo_search_new(protein, match, pa):
    aligner = Align.PairwiseAligner()
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    aligner.substitution_matrix = blosum62
    
    Lb = [0] * len(protein.terms)
    Ub = [match.frags_count[i] for i in match.frags_count]
    
    # generate individual by randomly select one from each TERM's match sequence
    def create_individual():
        '''
        individual = []
        for i in protein.terms:
            individual.append(round(random() * match.frags_count[i]))
        '''
        individual = [None] * len(protein.terms)
        for i in range(len(protein.terms)):
            nb = np.random.rand()
            if nb <= 0.75:
                individual[i] = round(random() * match.frags_count[protein.terms[i]])
            else:
                individual[i] = match.frags_count[protein.terms[i]]

        return individual
    
    # create initial population
    # input: size of population
    # output: a 2-dimensional numpy ndarray, each row is an individual
    def initial_population(popSize):
        population = np.zeros(shape = (popSize,len(protein.terms)),dtype = int)   
        for i in range(popSize):
           population[i] = create_individual()       
        return population
    
    # calculate score for each individual based on the amino acids alignment
    # input: individual
    # output: total alignment score of the input individual
    def compare_aa(individual):   
        frag_num = dict(zip(protein.terms, individual))
        frag_seq = {}
        for i in protein.terms:
                frag_seq[i] = match.select_frag(i, frag_num[i])      

        score = 0  
        for edge in protein.graph.edges:
            if type(frag_seq[edge[0]]) == np.ndarray and type(frag_seq[edge[1]]) == np.ndarray:
                 for pos in protein.graph.edges[edge]['sameAA']:
                    u_aa = frag_seq[edge[0]][pos[0]]
                    v_aa = frag_seq[edge[1]][pos[1]]
                    score += aligner.score(u_aa, v_aa)
                    
        return score/len(protein.terms)

    def term_count(individual):
        sel_frag = dict(zip(protein.terms, individual))    
        is_null = [sel_frag[i] < match.frags_count[i] for i in protein.terms]
        return np.sum(is_null)/len(protein.terms)
    
    def energy(individual):
        return compare_aa(individual) - 10 * term_count(individual) 
    
    # Simple bounds of the search domain
    def simple_bounds(s, Lb, Ub):
        for i in range(len(s)):
            if s[i] < Lb[i]:
                s[i] = Lb[i]
            if s[i] > Ub[i]:
                s[i] = Ub[i]
        s = s.astype(np.int)
        return s
    
    # Get cuckoos by random walk
    # Levy flights by Mantegna's algorithm
    # the difference factor (s - best) means that when the solution is the best solution, 
    # it remains unchanged
    def get_cuckoos(nest, best, Lb, Ub):
        beta = 1.5
        sigma_u = (math.gamma(1 + beta) * math.sin(math.pi * beta / 2) / (
                    math.gamma((1 + beta) / 2) * beta * (2 ** ((beta - 1) / 2)))) ** (1 / beta)
        sigma_v = 1
        for i in range(len(nest)):
            s = nest[i]
            u = np.random.normal(0, sigma_u, 1)
            v = np.random.normal(0, sigma_v, 1)
            step = u / ((abs(v)) ** (1 / beta))            
            stepsize = 0.01 * step * (s - best)   
            s = s + stepsize * np.random.randn(len(s))  
            nest[i] = simple_bounds(s, Lb, Ub)
        return nest
    
    # Find the current best nest
    def get_best_nest(nest, newnest, fitness):
        for j in range(len(nest)):
            fnew = energy(newnest[j])
            if fnew >= fitness[j]:
                fitness[j] = fnew
                nest[j] = newnest[j]        
        fmax = max(fitness.values())     
        bestnest =  max(fitness, key = fitness.get)
        return fmax, bestnest, nest, fitness
    
    # Replace some nests by constructing new solutions/nests
    # A fraction of worse nests are discovered with a probability pa
    def empty_nests(nest, Lb, Ub, pa):
        K = np.random.uniform(0, 1, size = nest.shape) > pa
        nest1 = nest.copy()
        nest2 = nest.copy()
        np.random.shuffle(nest1)
        np.random.shuffle(nest2)
        stepsize = np.random.uniform(0,1) * (nest1 - nest2)
        new_nest = nest + (stepsize * K).astype(int)
        for i in range(len(new_nest)):
            new_nest[i] = simple_bounds(new_nest[i], Lb, Ub)
            
        return new_nest

    # restore the sequence to letter form
    def restore_seq(individual):
        possible_aa = list([] for i in protein.terms)            
        candidate = dict(zip(protein.terms, possible_aa)) # inverse tells us for each residue, which protein.terms include it
        frag_num = dict(zip(protein.terms, individual)) # predict represents the choice of fragment for each TERM
        frag_seq = {}
        for i in protein.terms:
            frag_seq[i] = match.select_frag(i, frag_num[i]) 
    
        for pos in protein.inverse:
            for term in protein.inverse[pos]:
                if type(frag_seq[term]) == np.ndarray:
                    indice = protein.neighbors[term].index(pos)
                    candidate[pos].append(frag_seq[term][indice])

        possible_seq = ''
        for i in protein.terms:
            if candidate[i] != []:
                possible_seq += Counter(candidate[i]).most_common(1)[0][0]
            else:
                possible_seq += '-'
                continue
    
        return possible_seq
    
    # extract sequence from PDB
    def real_seq():
        structure = PDBParser().get_structure(protein.protein_id, protein.protein_id + '.pdb')
        ppb=PPBuilder() 
        seq = ''
        for pp in ppb.build_peptides(structure):
            seq += pp.get_sequence()
        return seq
    original = real_seq()
    
    # Random initial solutions
    nest = initial_population(100)
    # Get the current best
    fitness = {}
    for i in range(len(nest)):
        fitness[i] = energy(nest[i])
    fmax, bestnest, nest, fitness = get_best_nest(nest, nest, fitness)
    
    fitscore = []
    alignscore = []
    # Starting iterations
    for i in range(50):
        # Generate new solutions (but keep the current best)
        nest_a = nest.copy()
        new_nest = get_cuckoos(nest_a, bestnest, Lb, Ub)
        fnew, bestnest, nest, fitness = get_best_nest(nest, new_nest, fitness)
        # Discovery and randomization
        nest_b = nest.copy()
        new_nest = empty_nests(nest_b, Lb, Ub, pa)
    
        # Evaluate this solution
        fnew, bestnest, nest, fitness = get_best_nest(nest, new_nest, fitness)
        fitscore.append(fnew)
        predict = restore_seq(nest[bestnest])
        alignscore.append(aligner.score(predict, original))
    
    plt.plot(fitscore)
    plt.savefig(protein.protein_id + '_fitness.jpg')
    plt.close()
    plt.plot(alignscore)
    plt.savefig(protein.protein_id + '_align.jpg')
    plt.close()
    #return(score)

