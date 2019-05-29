import time
from function import *
from seq_pred import *

def runfunc(func):
    start = time.time()
    func
    end = time.time()
    print(end,start)
    #return end-start

fragments = node_attributes(match_sequence)


runfunc(createIndividual(fragments))
runfunc(initialPopulation(100, fragments))
runfunc(getFitness(pop[0], fragments, G))
runfunc(selection(pop, 50, fragments, G))
runfunc(crossover(pop[0], pop[1], 3))
runfunc(crossoverPopulation(pop, 3))
runfunc(mutate(pop[0], 0.03, fragments))
runfunc(mutatePopulation(pop, 0.03, fragments))
runfunc(getFitness(pop[0], fragments, G))
runfunc(nextGeneration(pop, 50, 3, 0.03, fragments, G))
runfunc(geneticAlgorithm(100, 50, 3, 0.03, fragments, 5, G))

