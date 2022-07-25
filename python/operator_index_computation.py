# -*- coding: utf-8 -*-
"""
Created on Sat Jul 23 00:50:57 2022

@author: rosha
"""

import csv
import numpy as np
from pygmo import hypervolume

### Parameters
assigning_problem = False
random_mode = 1 # 1 - only random data, 2 - only epsilon MOEA, 3 - both

num_runs = 10

num_archs = 300

save_path = "C:\\Users\\rosha\\Documents\\SEAK Lab Github\\VASSAR\\VASSAR_exec_heur\\results" # for laptop
#save_path = "C:\\SEAK Lab Github\\VASSAR\\VASSAR_exec_heur\\results" # for workstation

filepath_prob = "Assigning\\"
if (not assigning_problem):
    filepath_prob = "Partitioning\\"
    
filepath_mode = "Random\\"
if assigning_problem:
    filename_mode = "random_assign_operator_index_"
else:
    filename_mode = "random_partition_operator_index_"
if (random_mode == 2):
    filepath_mode = "Epsilon MOEA\\"
    if assigning_problem:
        filename_mode = "EpsilonMOEA_emoea_ClimateCentric_assign_operator_index_"
    else:
        filename_mode = "EpsilonMOEA_emoea_ClimateCentric_partition_operator_index_"
elif (random_mode == 3):
    filepath_moea = "Epsilon MOEA\\"
    if assigning_problem:
        filename_moea = "EpsilonMOEA_emoea_ClimateCentric_assign_operator_index_"
    else:
        filename_moea = "EpsilonMOEA_emoea_ClimateCentric_partition_operator_index_"
        
### Useful functions
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

def read_csv_run(filepath, filename, run_num, run_mode, n_archs):
    full_filename = filepath + filename + run_num + '.csv'
    
    with open(full_filename,newline='') as csvfile:
        data = [row for row in csv.reader(csvfile)]
        
        if (run_mode == 2) or (run_mode == 3):
            nfes = np.zeros(len(data)-1)
        
        science = np.zeros(len(data)-1)
        cost = np.zeros(len(data)-1)
        
        science_instrdc = np.zeros(len(data)-1)
        cost_instrdc = np.zeros(len(data)-1)
    
        science_instrorb = np.zeros(len(data)-1)
        cost_instrorb = np.zeros(len(data)-1)
    
        science_interinstr = np.zeros(len(data)-1)
        cost_interinstr = np.zeros(len(data)-1)
        
        science_packeff = np.zeros(len(data)-1)
        cost_packeff = np.zeros(len(data)-1)
        
        science_spmass = np.zeros(len(data)-1)
        cost_spmass = np.zeros(len(data)-1)
        
        science_instrsyn = np.zeros(len(data)-1)
        cost_instrsyn = np.zeros(len(data)-1)
        
        buffer = 0
        if (run_mode == 2) or (run_mode == 3):
            buffer = 1
        
        valid_count = 0
        for i in range(len(data)-1):
            
            if (run_mode == 2) or (run_mode == 3):
                nfe = list(map(float,data[i+1][buffer]))
            
            data_float = list(map(float,data[i+1][buffer+1:buffer+2]))
            data_float_instrdc = list(map(float,data[i+1][buffer+4:buffer+5]))
            data_float_instrorb = list(map(float,data[i+1][buffer+7:buffer+8]))
            data_float_interinstr = list(map(float,data[i+1][buffer+9:buffer+10]))
            data_float_packeff = list(map(float,data[i+1][buffer+11:buffer+12]))
            data_float_spmass = list(map(float,data[i+1][buffer+13:buffer+14]))
            data_float_instrsyn = list(map(float,data[i+1][buffer+15:buffer+16]))
            
            if (run_mode == 2) or (run_mode == 3):
                data_all = nfes + data_float + data_float_instrdc + data_float_instrorb + data_float_interinstr + data_float_packeff + data_float_spmass + data_float_instrsyn
            else:
                data_all = data_float + data_float_instrdc + data_float_instrorb + data_float_interinstr + data_float_packeff + data_float_spmass + data_float_instrsyn
            
            if (any(np.isnan(np.array(data_all))) or any(np.isinf(np.array(data_all)))):
                continue
            
            if (run_mode == 2) or (run_mode == 3):
                nfes[valid_count] = nfe
            
            science[valid_count] = data_float[0]
            cost[valid_count] = data_float[1]
            
            science_instrdc[valid_count] = data_float_instrdc[0]
            cost_instrdc[valid_count] = data_float_instrdc[1]
            
            science_instrorb[valid_count] = data_float_instrorb[0]
            cost_instrorb[valid_count] = data_float_instrorb[1]
            
            science_interinstr[valid_count] = data_float_interinstr[0]
            cost_interinstr[valid_count] = data_float_interinstr[1]
            
            science_packeff[valid_count] = data_float_packeff[0]
            cost_packeff[valid_count] = data_float_packeff[1]
            
            science_spmass[valid_count] = data_float_spmass[0]
            cost_spmass[valid_count] = data_float_spmass[1]
            
            science_instrsyn[valid_count] = data_float_instrsyn[0]
            cost_instrsyn[valid_count] = data_float_instrsyn[1]
            
            
            
            