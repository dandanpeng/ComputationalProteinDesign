import re, copy
import numpy as np
import operator 
from Bio import Align
from Bio.SubsMat.MatrixInfo import blosum62
from random import sample, choice, random
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

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
# input: fragments is a dictonary (keys: TERMs' names, values: matched sequence) 
# output:  list, each element represents the index of randomly chosen sequence
def createIndividual(fragments):
    individual = []
    for i in fragments.seq:
        individual.append(round(random() * len(fragments.seq[i])))
    return individual

# create initial population, 
# input: size of population, fragments
# output: a 2-dimensional numpy array, each row is an individual
def initialPopulation(popSize, fragments):
    population = np.zeros(shape = (popSize,len(fragments.seq)),dtype = int)
    
    for i in range(popSize):
       population[i] = createIndividual(fragments)
    return population

# select elites from children
# input: population, size of elites, fragmants, graph that represents topolpgy pf protein
# output: individuals who have high score
def selection(population, eliteSize, fragments, G):
    fitnessResults = {}
    for i in range(len(population)):
        fitnessResults[i] = get_fitness(population[i],fragments, G)
    sortedResults = sorted(fitnessResults.items(), key = operator.itemgetter(1), reverse = True)
    elites = [i[0] for i in sortedResults[:eliteSize]]   
    non_elites =  list(np.random.choice(range(len(population)), size = eliteSize))
    matingpool = np.append(population[elites,:],population[non_elites,:],axis = 0)
    return matingpool     

# simulate crossover process between two parents
# input: two parents (list), and num_points (the number of points between which crossover happens)
# output: child produced by parent1 and parent2       
def crossover(parent1, parent2,num_points):
    points = sample([i for i in range(len(parent1))], num_points)
    points = list(set(points) | {0,len(parent1)})
    points.sort()
    child = []
    i = 0  
    names = locals()
    while i < len(points)-1: 
        names['child' + str(i)] = choice([parent1[points[i]:points[i+1]],parent2[points[i]:points[i+1]]])
        child = child + list(names.get('child'+str(i)))
        i+=1  
    return child

# simulate crossover process among population, randomly select two individuals as parents each time
# input: matingpool (candidate parents), popSize (size of children), num_points
# output: children (2-dimensional array, each row is an individual)
def crossoverPopulation(matingpool,num_points):
    children = np.zeros(shape = (matingpool.shape[0], matingpool.shape[1]), dtype = int)
    for i in range(len(children)):
        parents = sample(range(len(children)),2)
        children[i] = crossover(matingpool[parents[0]],matingpool[parents[1]], num_points)
    return children

# simulate gene mutation on individual
# input: individual, mutationRate (rate of mutation), fragments
# output: individual (mutated individual)        
def mutate(individual,mutationRate,fragments):
    mutationNumber = int(mutationRate * len(individual))
    mutationPosition = sample(range(len(individual)),mutationNumber)
    for i in mutationPosition:
        individual[i] = choice(range(len(fragments.seq['A'+str(i+1)])))
    return individual

# simulate gene mutation on population
# input: population, mutationRate, fragments
# output: mutatePop (mutated population)
def mutatePopulation(population, mutationRate, fragments):
    mutatedPop = copy.deepcopy(population)
    for ind in range(len(population)):
        mutatedPop[ind] = mutate(mutatedPop[ind],mutationRate,fragments)
        return mutatedPop

# calculate score for each individual based on the amino acids    
# input: individual, fragments, G (graph)
# output: score of the input individual
def get_fitness(individual,fragments, G):
    aligner = Align.PairwiseAligner()
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    aligner.substitution_matrix = blosum62
    
    keys = sort_string(fragments.seq.keys())
    rand_seq = dict(zip(keys,individual))
    score = 0
    
    for i in G.edges:
        if rand_seq[i[0]] < len(fragments.seq[i[0]]) and rand_seq[i[1]] < len(fragments.seq[i[1]]):
            for j in G.edges[i]['sameAA']:
                u_seq = (fragments.select(i[0],rand_seq[i[0]]))[j[0]]
                v_seq = (fragments.select(i[1],rand_seq[i[1]]))[j[1]]
                score += aligner.score(u_seq,v_seq)        
    return score



def nextGeneration(population, eliteSize, num_points, mutationRate, fragments, G):
    matingpool = selection(population, eliteSize, fragments, G)
    children = crossoverPopulation(matingpool, num_points)
    nextGeneration = mutatePopulation(children, mutationRate, fragments)
    return nextGeneration



def avScorePop(population, fragments, G):
    score = []
    for i in range(population.shape[0]):
        score.append(get_fitness(population[i],fragments,G))
    return np.mean(score)



def geneticAlgorithm(popSize, eliteSize, num_points, mutationRate, fragments, generations,G):
    pop = initialPopulation(popSize, fragments)
    print("Initial population finished")
    
    for i in range(generations):
        pop = nextGeneration(pop, eliteSize, num_points, mutationRate, fragments, G)
      
    return pop



def geneticAlgorithmPlot(popSize, eliteSize, num_points, mutationRate, fragments, generations,G):
    pop = initialPopulation(popSize, fragments)
    
    avg = []
    for i in range(generations):
        pop = nextGeneration(pop, eliteSize, num_points, mutationRate, fragments, G)
        avg.append(avScorePop(pop,fragments,G))
    
    plt.plot(avg)
    plt.ylabel('Average score of each generation')
    plt.xlabel('Generation')
    plt.savefig("/cluster/home/pengd/project/GAplot.jpg")  


def restore_seq(keys, individual, neighbors,inverse):
    possible_res = list([] for i in keys)            
    candidate = dict(zip(keys,possible_res)) # inverse tells us for each residue, which TERMs include it
    predict = dict(zip(keys,list(individual))) # predict represents the choice of fragment for each TERM
    
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
    
    return possible_seq