import operator, argparse
import numpy as np

#import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as plt

from seqpred import *
from random import sample, choice, random
from collections import Counter
from Bio import Align
from Bio.PDB.PDBParser import PDBParser
from Bio.SubsMat.MatrixInfo import blosum62
from Bio.PDB.Polypeptide import PPBuilder

par = argparse.ArgumentParser()
par.add_argument('--p', required = True, help = 'protein ID')
args = par.parse_args()

protein_id = args.p

alpha = 100
neighbors, keys, inverse, overlap, match, frag, frags_count, G = get_information(protein_id)
aligner = Align.PairwiseAligner()
aligner.open_gap_score = -10
aligner.extend_gap_score = -0.5
aligner.substitution_matrix = blosum62    

######### Genetic Algorithm #############
# input: list frag_count, which stores the number of fragments each TERM has
# output:  list, each element represents the index of randomly chosen sequence
def create_individual():
    individual = []
    for i in keys:
        individual.append(round(random() * frags_count[i]))   
    return individual

# create initial population, 
# input: size of population, fragments
# output: a 2-dimensional numpy array, each row is an individual
def initial_population(popSize):
    population = np.zeros(shape = (popSize,len(keys)),dtype = int)   
    for i in range(popSize):
       population[i] = create_individual()       
    return population


# calculate score for each individual based on the amino acids    
# input: individual, fragments, G (graph)
# output: score of the input individual

def compare_aa(individual):   
    sel_frag = dict(zip(keys, individual))
    sel_frag_seq = {}
    for i in keys:
        if sel_frag[i] < frags_count[i]:
            sel_frag_seq[i] = frag.select(i, sel_frag[i])      
    
    score = 0  
    for edge in G.edges:
        if sel_frag[edge[0]] < frags_count[edge[0]] and sel_frag[edge[1]] < frags_count[edge[1]]:
             for pos in G.edges[edge]['sameAA']:
                u_aa = sel_frag_seq[edge[0]][pos[0]]
                v_aa = sel_frag_seq[edge[1]][pos[1]]
                score += aligner.score(u_aa, v_aa)
    return score/len(keys)


def term_count(individual):
    sel_frag = dict(zip(keys, individual))    
    is_null = [sel_frag[i] < frags_count[i] for i in keys]
    return np.sum(is_null)/len(keys)

def energy(individual):
    return compare_aa(individual) - alpha * term_count(individual) 
            
# select elites from children
# input: population, size of elites, fragmants, graph that represents topolpgy pf protein
# output: individuals who have high score
def selection(population, eliteSize):
    fitnessResults = {}
    for i in range(len(population)):
        fitnessResults[i] = energy(population[i])

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
def mutate(individual, mutationRate):
    mutationNumber = int(mutationRate * len(individual))
    mutationPosition = sample(range(len(individual)),mutationNumber)
    for i in mutationPosition:
        individual[i] = choice(range(frags_count[keys[i]]))
    return individual

# simulate gene mutation on population
# input: population, mutationRate, fragments`
# output: mutatePop (mutated population)
def mutate_population(population, mutationRate):
    mutatedPop = population.copy()
    for ind in range(len(population)):
        mutatedPop[ind] = mutate(mutatedPop[ind], mutationRate)
        return mutatedPop
        

def next_generation(population, eliteSize, num_points, mutationRate):
    matingpool = selection(population, eliteSize)
    children = crossover_population(matingpool, num_points)
    nextGeneration = mutate_population(children, mutationRate)
    return nextGeneration



def genetic_algorithm(popSize, eliteSize, num_points, mutationRate, generations):
    pop = initial_population(popSize)
    
    for i in range(generations):
        pop = next_generation(pop, eliteSize, num_points, mutationRate)
      
    return pop

# restore to letter form
def restore_seq(individual):
    possible_res = list([] for i in keys)            
    candidate = dict(zip(keys, possible_res)) # inverse tells us for each residue, which TERMs include it
    sel_frag = dict(zip(keys, individual)) # predict represents the choice of fragment for each TERM
    
    if_null = [sel_frag[i] >= frag.count(i) for i in keys]
    if_null_dic = dict(zip(keys, if_null))
    
    for pos in inverse:
        for term in inverse[pos]:
            if if_null_dic[term] == False:
                indice = neighbors[term].index(pos)
                candidate[pos].append(frag.select(term, sel_frag[term])[indice])

    possible_seq = ''
    for i in keys:
        if candidate[i] != []:
            possible_seq += Counter(candidate[i]).most_common(1)[0][0]
        else:
            possible_seq += '-'
            continue
    
    return possible_seq

# extract sequence from PDB
def real_seq():
    structure = PDBParser().get_structure(protein_id, protein_id + '.pdb')
    ppb=PPBuilder() 
    seq = ''
    for pp in ppb.build_peptides(structure):
        seq += pp.get_sequence()
    return seq

pop = genetic_algorithm(100, 50, 2, 0.03, 10)
pred_seq = restore_seq(pop[0])
original_seq = real_seq()

print(aligner.score(original_seq, pred_seq))

############# Plot function ###############

def len_frag(population):
    len_frags = [frag.seq[i].shape[1] for i in sort_string(frag.seq.keys())]
    nb_frags = [frag.seq[i].shape[0] for i in sort_string(frag.seq.keys())]
    
    len_frag = 0
    for individual in population:
        diff = individual - nb_frags
        for i in range(len(diff)):
            if diff[i] < 0:
                len_frag += len_frags[i]
                
    return len_frag/(len(population) * len(frag.seq))
    
    
    

def plot(popSize, eliteSize, num_points, mutationRate, generations):
    pop = initial_population(popSize)
    
    fitness = []
    score = []
    length = []
    
    for i in range(generations): 
        pop = next_generation(pop, eliteSize, num_points, mutationRate)        
        #nb_pop_frag = term_count(pop[0:10], keys, frags_count)
        #len_pop_frag = len_frag(pop[0:10],frag)
        fitness.append(energy(pop[0]))
        pred_seq = restore_seq(pop[0])
        score.append(aligner.score(original_seq, pred_seq))
        #number.append(term_count(pop[0]))
        #length.append(len_pop_frag)
        
    plt.plot(fitness)
    plt.ylabel('Energy of each generation\'s first elite')
    plt.xlabel('Generation')
    plt.savefig('/Users/pengdandan/Desktop/lab_rotation/LabRotation2/code/' + protein_id + '_energy.jpg")  
    plt.close()
    
    plt.plot(score)
    plt.ylabel('Align score of each generation\'s first elite')
    plt.xlabel('Generation')
    plt.savefig('/Users/pengdandan/Desktop/lab_rotation/LabRotation2/code/' + protein_id + 'alignScore.jpg")  
    plt.close()
    
    '''
    plt.plot(number)
    plt.ylabel('Average frag number of each generation')
    plt.xlabel('Generation')
    #plt.savefig("/cluster/home/pengd/project/test/nbplot.jpg")  
    plt.close()
    
    plt.plot(length)
    plt.ylabel('Average frag length of each generation')
    plt.xlabel('Generation')
    #plt.savefig("/cluster/home/pengd/project/test/lenplot.jpg")  
'''
#############       
        

