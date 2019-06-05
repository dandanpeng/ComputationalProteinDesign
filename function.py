import re, copy
import numpy as np
import operator 
import multiprocessing

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

from Bio import Align
from Bio.SubsMat.MatrixInfo import blosum62
from random import sample, choice, random
from collections import Counter



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


######### Genetic Algorithm #############
# input: list frag_count, which stores the number of fragments each TERM has
# output:  list, each element represents the index of randomly chosen sequence
def create_individual(frag):
    individual = []
    for i in sort_string(frag.seq.keys()):
        individual.append(round(random() * frag.count(i)))    
    return individual

# create initial population, 
# input: size of population, fragments
# output: a 2-dimensional numpy array, each row is an individual
def initial_population(popSize, frag):
    population = np.zeros(shape = (popSize,len(frag.seq)),dtype = int)
    
    for i in range(popSize):
       population[i] = create_individual(frag)
       
    return population


# calculate score for each individual based on the amino acids    
# input: individual, fragments, G (graph)
# output: score of the input individual

def compare_aa(individual, keys, frag, G):    
    aligner = Align.PairwiseAligner()
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    aligner.substitution_matrix = blosum62
    
    sel_frag = dict(zip(keys, individual))
    score = 0  
    
    for edge in G.edges:
        if sel_frag[edge[0]] < frag.count(edge[0]) and sel_frag[edge[1]] < frag.count(edge[1]):
            for pos in G.edges[edge]['sameAA']:
                u_aa = frag.select(edge[0],sel_frag[edge[0]])[pos[0]]
                v_aa = frag.select(edge[1],sel_frag[edge[1]])[pos[1]]
                score += aligner.score(u_aa, v_aa)

    return score

## parallel programming
def compare_aa(indivual, keys, frag, G):
        aligner = Align.PairwiseAligner()
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    aligner.substitution_matrix = blosum62
    
    sel_frag = dict(zip(keys, individual))
    score = 0  
        
    for i in sel_frag:
        if sel_frag[i] >= frag.count(i):
            len_term = frag.seq[i].shape[1]
            frag.seq[i] = np.append(frag.seq[i], [['-'] * len_term], axis = 0)
    
    edges1 = G.edges[0:331]
    edges2 = G.edges[331: 761]
    pos1 = [G.edges[i]['sameAA'] for i in edges1]
    pos2 = [G.edges[i]['sameAA'] for i in egdes2]
    
    pool = multiprocessing.Pool()
    res = pool.map(score, paramlist)

# restore to letter form
def restore_seq(keys, individual, neighbors, frag, G):
    possible_res = list([] for i in keys)            
    candidate = dict(zip(keys,possible_res)) # inverse tells us for each residue, which TERMs include it
    sel_frag = dict(zip(keys,individual)) # predict represents the choice of fragment for each TERM
           
    for i in sel_frag:
        if sel_frag[i] >= frag.count(i):
            len_term = frag.seq[i].shape[1]
            frag.seq[i] = np.append(frag.seq[i], [['-'] * len_term], axis = 0)
    
    for edge in G.edges:
        u_frag = frag.select(edge[0],sel_frag[edge[0]])
        v_frag = frag.select(edge[1],sel_frag[edge[1]])
        for pos in G.edges[edge]['sameAA']:
            aa = (neighbors[edge[0]])[pos[0]] ## aa is every amino acid shared by node u and v
            candidate[aa].append(u_frag[pos[0]])
            candidate[aa].append(v_frag[pos[1]])

    possible_seq = ''
    for i in keys:
        if candidate[i] != []:
            possible_seq += Counter(candidate[i]).most_common(1)[0][0]
        else:
            possible_seq += '-'
            continue
    
    return possible_seq


def individual_gap(keys, individual, neighbors, frag, G):
    seq = restore_seq(keys, individual, neighbors, frag, G)
    return seq.count('-')


def nb_frag(individual, frag):
    nb_frags = [frag.count[i]) for i in sort_string(frag.seq.keys())]
    diff = individual - nb_frags
    return (diff < 0).sum(0)

def energy_func(individual, keys, frag, G, neighbors):
    return compare_aa(individual, keys, frag, G) 
            - individual_gaps(keys, individual, neighbor, frag, G)
            - nb_frag(individual, frag) 
            
# select elites from children
# input: population, size of elites, fragmants, graph that represents topolpgy pf protein
# output: individuals who have high score
def selection(population, eliteSize, frag, G, keys):
    fitnessResults = {}
    
    for i in range(len(population)):
        fitnessResults[i] = get_fitness(population[i], keys, frag, G)

    sortedResults = sorted(fitnessResults.items(), key = operator.itemgetter(1), reverse = True)
    
    
    elites = [i[0] for i in sortedResults[:eliteSize]]   
    non_elites =  list(np.random.choice(range(len(population)), size = eliteSize))
    matingpool = np.append(population[elites,:],population[non_elites,:],axis = 0)
    return matingpool     


# simulate crossover process between two parents
# input: two parents (list), and num_points (the number of points between which crossover happens)
# output: child produced by parent1 and parent2       
def crossover(parent1, parent2, num_points):
    points = sample(range(len(parent1)), num_points)
    points = list(set(points) | {0, len(parent1)})
    points.sort()
    child = []
    for i in range(len(points)-1):
        provider = choice([parent1, parent2])
        child.extend(provider[points[i]:points[i+1]])
    
    return child
    
    
# simulate crossover process among population, randomly select two individuals as parents each time
# input: matingpool (candidate parents), popSize (size of children), num_points
# output: children (2-dimensional array, each row is an individual)
def crossover_population(matingpool, num_points):
    children = np.zeros(shape = (matingpool.shape[0], matingpool.shape[1]), dtype = int)
    for i in range(len(children)):
        parents = sample(list(matingpool),2)
        children[i] = crossover(parents[0], parents[1], num_points)
    return children

# simulate gene mutation on individual
# input: individual, mutationRate (rate of mutation), fragments
# output: individual (mutated individual)        
def mutate(individual, mutationRate, frag, keys):
    mutationNumber = int(mutationRate * len(individual))
    mutationPosition = sample(range(len(individual)),mutationNumber)
    for i in mutationPosition:
        individual[i] = choice(range(frag.count(keys[i])))
    return individual

# simulate gene mutation on population
# input: population, mutationRate, fragments`
# output: mutatePop (mutated population)
def mutate_population(population, mutationRate, frag, keys):
    mutatedPop = population.copy()
    for ind in range(len(population)):
        mutatedPop[ind] = mutate(mutatedPop[ind], mutationRate, frag, keys)
        return mutatedPop
        

def next_generation(population, eliteSize, num_points, mutationRate, frag, G, keys):
    matingpool = selection(population, eliteSize, frag, G, keys)
    children = crossover_population(matingpool, num_points)
    nextGeneration = mutate_population(children, mutationRate, frag, keys)
    return nextGeneration



def genetic_algorithm(popSize, eliteSize, num_points, mutationRate, frag, generations, G, keys):
    pop = initial_population(popSize, frag, keys)
    
    for i in range(generations):
        pop = next_generation(pop, eliteSize, num_points, mutationRate, frag, G, keys)
      
    return pop


############# Plot function ###############
def av_sc(population, frag, G):
    score = []
    for i in range(population.shape[0]):
        score.append(getFitness(population[i],frag,G))
    return np.mean(score)



def nb_frag(population, frag):
    nb_frags = [len(frag.seq[i]) for i in sort_string(frag.seq.keys())]
    nb_frag = 0
    
    for individual in population:
        diff = individual - nb_frags
        nb_frag += (diff < 0).sum(0)
    
    return nb_frag/len(population)



def len_frag(population, frag):
    len_frags = [frag.seq[i].shape[1] for i in sort_string(frag.seq.keys())]
    nb_frags = [frag.seq[i].shape[0] for i in sort_string(frag.seq.keys())]
    
    len_frag = 0
    for individual in population:
        diff = individual - nb_frags
        for i in range(len(diff)):
            if diff[i] < 0:
                len_frag += len_frags[i]
                
    return len_frag/(len(population) * len(frag.seq))
    


def plot(popSize, eliteSize, num_points, mutationRate, frag, generations, G):
    pop = initial_population(popSize, frag, keys)
    
    score = []
    number = []
    length = []
    for i in range(generations): 
        pop = next_generation(pop, eliteSize, num_points, mutationRate, frag, G)        
        nb_pop_frag = nb_frag(pop[0:10],frag)
        len_pop_frag = len_frag(pop[0:10],frag)
        score.append(av_sc(pop[0:10],frag,G))
        number.append(nb_pop_frag)
        length.append(len_pop_frag)
        
    plt.plot(score)
    plt.ylabel('Average score of each generation')
    plt.xlabel('Generation')
    #plt.savefig("/cluster/home/pengd/project/test/scplot.jpg")  
    plt.close()
    
    plt.plot(number)
    plt.ylabel('Average frag number of each generation')
    plt.xlabel('Generation')
    #plt.savefig("/cluster/home/pengd/project/test/nbplot.jpg")  
    plt.close()
    
    plt.plot(length)
    plt.ylabel('Average frag length of each generation')
    plt.xlabel('Generation')
    #plt.savefig("/cluster/home/pengd/project/test/lenplot.jpg")  

        
        

