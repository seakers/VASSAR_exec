# -*- coding: utf-8 -*-
"""
Hypervolume function testing

@author: rosha
"""
from pygmo import hypervolume
import numpy as np
import matplotlib.pyplot as plt

def compute_pareto_front(population):
    pop_size = len(population)
    obj_num = 2

    domination_counter = [0] * pop_size

    for i in range(pop_size):
        for j in range(i+1, pop_size):
            # check each objective for dominance
            dominate = [0] * obj_num
            for k in range(obj_num):
                if population[i][k] > population[j][k]:
                    dominate[k] = 1
                elif population[i][k] < population[j][k]:
                    dominate[k] = -1
            if -1 not in dominate and 1 in dominate:
                domination_counter[i] += 1
            elif -1 in dominate and 1 not in dominate:
                domination_counter[j] += 1

    pareto_solutions = []
    for i in range(len(domination_counter)):
        if domination_counter[i] == 0:
            pareto_solutions.append(population[i])
    return pareto_solutions

def compute_hv(population):
    array_archs = np.zeros((len(population), 2))
    for i in range(len(population)):
        array_archs[i] = population[i]
    hv_object = hypervolume(array_archs)
    hv = hv_object.compute([1.1,1.1])/1.1**2
    return hv

def plot_test_population(population):
    plt.figure()
    x_pop = [x[0] for x in population]
    y_pop = [x[1] for x in population]
    plt.scatter(x_pop, y_pop)
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.show()

#test_population = [[0, 1.1], [1.1, 0], [0.55, 0.55]]
test_population = [[0.9, 0.5], [0.6, 0.6], [0.1, 0.4], [0.8, 0.2], [0.5, 0.5], [0.3, 0.6]]
#test_population = [[1, 1]] # worst point
#test_population = [[0, 0]] # best point

plot_test_population(test_population)

pareto_front = compute_pareto_front(test_population)
print('Pareto Front = ', pareto_front)

hv_val = compute_hv(pareto_front)
print('Hypervolume for the test population = ', hv_val)
