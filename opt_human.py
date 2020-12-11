from deap import base, creator, tools
import file_comparison
import json
import copy
import itertools as it
import os
import numpy as np
from rule_creator import select_rules_human_guided
import ast
from time import time
from joblib import Parallel, delayed
import multiprocessing
import random


#load reference data
fn_original = 'EMT_paper.json'
symbols = ['NICD', 'Notch', 'TP53', 'TP63TP73', 'miRNA', 'EMTreg', 'ECM', 'DNAdam']


(N, idx_freq) = file_comparison.map_freqs(fn_original)
N = 6

bounds = np.array([(0, 100)]*N)
x0 = np.random.rand(N)*100

pop_size = 150
generation = 0
counter = 0
parallel_sims = 1


def randomized_runs(perc):
    str_rules, rule_list, file_name = select_rules_human_guided(perc, symbols)
    exe_file = '{}.bin'.format(file_name)
    json_file = '{}.json'.format(file_name)
    os.system('time g++ -DOUTPUT_FILE=\\"{2}\\" /home/martina/jsoncpp-master/jsoncpp-master/dist/jsoncpp.cpp -I/home/martina/jsoncpp-master/jsoncpp-master/dist/json  -O3 -fopenmp -x c++ {0} -o {1} && {1}'.format(file_name,exe_file,json_file))

    with open(json_file) as fs:
        freq_multiple_ss = json.load(fs)
    fs.close()

    os.remove(file_name)
    os.remove(exe_file)
    os.remove(json_file)

    #print("file comparison: ", file_comparison.comparator(freq_multiple_ss, fn_original))
    #print("freq_multiple_ss: ", freq_multiple_ss)
    #print("str_rules: ", str_rules)
    #exit()

    return file_comparison.comparator(freq_multiple_ss, fn_original), freq_multiple_ss, str_rules, rule_list

def fitness(perc_to_change):
    global counter, generation

    #print('here')
    # take params and write file
    # ==========================
    with open(fn_original) as f2:
        orig_data = json.load(f2)

    number_of_combinations = 6


    # run SSA
    # =======
    rms = 100.0

    #num_cores = int(multiprocessing.cpu_count()/2)
    num_cores = 1
    print('using {} cores'.format(num_cores))
    a = time()
    #print("perc_to_change", perc_to_change)
    #exit()
    #results = Parallel(n_jobs=num_cores)(delayed(lambda: randomized_runs(perc_to_change))() for i in range(parallel_sims))
    rms_new, freq_multiple_ss, str_rules, rule_list = randomized_runs(perc_to_change)
    #print(results)
    #exit()
    print('parallel region took {}'.format(time()-a))
    #print(results)

    fs_fitrules = open("optimization_output/fitness-{}-{}.txt".format(str(generation).zfill(3),str(counter).zfill(3)),"w")

    #for rms_new, freq_multiple_ss, str_rules, rule_list in results:
    rms_file.write('  {}\n'.format(rms_new))
    rms_file.flush()

    fs_fitrules.write('{}\n'.format(rms_new))
    fs_fitrules.write('{}'.format(str_rules))
    for freq_iv in freq_multiple_ss:
        fs_fitrules.write('{} -> '.format(freq_iv[0]))
        is_first = True
        for ss_freq in freq_iv[1]:
            if is_first == True:
                fs_fitrules.write('{} {}\n'.format(ss_freq[0], ss_freq[1]))
                is_first = False
            else:
                fs_fitrules.write('{:>15} {} {}\n'.format('', ss_freq[0], ss_freq[1]))
    for idx, rule in zip(range(1,len(rule_list)+1),rule_list):
        fs_fitrules.write("{}: {}\n".format(idx,rule))
    fs_fitrules.write('\n\n')

    if rms_new < rms:
        rms = rms_new

    fs_fitrules.close()
    print(generation, counter)
    #os.rename("modified_freq.json", "optimization_output/run-{}-{}.json".format(str(generation).zfill(3), str(counter).zfill(3)))

    rms_file.write('{}\t{}\t{}\n'.format(counter,rms, perc_to_change))

    counter += 1
    if counter % pop_size == 0:
        counter = 0
        generation += 1
        rms_file.write('\n')

    rms_file.flush()

    return rms,


# using DEAP optimizer
creator.create("FitnessMin", base.Fitness, weights=(-1.0,))
creator.create("Individual", list, fitness=creator.FitnessMin)

toolbox = base.Toolbox()
#initializer for the starting generation - each individual starts with this setup
toolbox.register("attr_float", random.random)  # the frequency starts with a random number between 0 and 100
toolbox.register("individual", tools.initRepeat, creator.Individual, toolbox.attr_float, N) # an individual has a list of N random numbers that will be changed by the optimizer
toolbox.register("population", tools.initRepeat, list, toolbox.individual)

#set up genetic optimizer functions
toolbox.register("evaluate", fitness)
#toolbox.register("mate", tools.cxTwoPoint)
toolbox.register("mate", tools.cxUniform, indpb=0.05)
#toolbox.register("mutate", tools.mutFlipBit, indpb=0.05)
toolbox.register("mutate", tools.mutPolynomialBounded, low=0.0, up=1.0, eta=20.0, indpb=0.1)
#toolbox.register("mutate", tools.mutUniformInt, low=0, up=1, indpb=0.1)
toolbox.register("select", tools.selTournament, tournsize=3)


def do_the_opt_dance():
    pop = toolbox.population(n=pop_size)  # note: here the population is only between 0 and 1 - this will be multiplied by 100 at the assignment in the fitness fct
    print("pop: ", pop)
    fitnesses = list(map(toolbox.evaluate, pop))  # map evaluation fct with every individual
    #print("fitnesses: " , fitnesses)
    #exit()

    for ind, fit in zip(pop, fitnesses):
        ind.fitness.values = fit

    CXPB = 0.5  # probability with which two individuals are crossed
    MUTPB = 0.5  # probability for mutating and individual

    fits = [ind.fitness.values[0] for ind in pop]

    # Variable keeping track of the number of generations
    g = 0
    # Begin the evolution
    while min(fits) > 0 and g < 100:
        # A new generation
        g = g + 1
        print("-- Generation %i --" % g)

        # evolution: selecting, mating and mutating the individuals in the population
        # step 1: Select the next generation individuals
        offspring = toolbox.select(pop, len(pop))
        # Clone the selected individuals
        offspring = list(map(toolbox.clone, offspring)) # offspring list is an exact copy of the selected individuals. toolbox.clone ensures that we don't use a reference to the individuals but a completely independent instance

        # simple example, step 2&3: Apply crossover and mutation on the offspring
        for child1, child2 in zip(offspring[::2], offspring[1::2]):
            if random.random() < CXPB:
                toolbox.mate(child1, child2) # mating with Uniform distribution -> set probability to mate to 25%
                del child1.fitness.values
                del child2.fitness.values
        # the mutation of the prdocut children with a certain probability of CXPB and MUTPB.
        for mutant in offspring:
            if random.random() < MUTPB:
                toolbox.mutate(mutant)
                del mutant.fitness.values #the del statement will invalidate the fitness of the modified offspring

        # Content of some of our offspring changed during last step -> need to re-evaluate their fitnesses (just map those offspring with fitnesses where marked invalid)
        # Evaluate the individuals with an invalid fitness
        invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
        #print(invalid_ind)
        #exit()
        fitnesses = map(toolbox.evaluate, invalid_ind)
        for ind, fit in zip(invalid_ind, fitnesses):
            ind.fitness.values = fit
        # Replace the old population by the offspring
        pop[:] = offspring





rms_file = open("rms_fit.txt", "w+")
do_the_opt_dance()
