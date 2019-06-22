import time
from function import *
#from seq_pred import *

fragments = node_attributes(match_sequence)

start = time.time()
create_individual(keys, frags_count)
end = time.time()
print("create_indiviaudl: " + str(end - start))

start = time.time()
pop = initial_population(100, keys, frags_count)
end = time.time()
print("initial_population: " + str(end - start))

start = time.time()
compare_aa(pop[0], keys, fragments, frags_count, G)
end = time.time()
print("compare_aa: " + str(end - start))

start = time.time()
term_count(pop[0], keys, frags_count)
end = time.time()
print("term_count: " + str(end - start))

start = time.time()
energy(pop[0], keys, frag, frags_count, G)
end = time.time()
print("energy: " + str(end - start))


start = time.time()
selection(pop, 50, fragments, G, keys)
end = time.time()
print("selection: " + str(end - start))

start = time.time()
crossover(pop[0], pop[1], 3)
end = time.time()
print("crossover: " + str(end - start))

start = time.time()
crossover_population(pop, 3)
end = time.time()
print("crossover_population: " + str(end - start))

start = time.time()
mutate(pop[0], 0.03, fragments, keys)
end = time.time()
print("mutate: " + str(end - start))

start = time.time()
mutate_population(pop, 0.03, fragments, keys)
end = time.time()
print("mutate_population: " + str(end - start))

start = time.time()
next_generation(pop, 50, 3, 0.03, fragments, G, keys)
end = time.time()
print("next_generation: " + str(end - start))

start = time.time()
genetic_algorithm(100, 50, 3, 0.03, fragments, frags_count, 1, G, keys)
end = time.time()
print("genetic_algorithm: " + str(end - start))

