# -*- coding: utf-8 -*-
"""
Created on Tue Jul 26 19:32:05 2022

@author: rosha
"""
import csv
import numpy as np
from pygmo import hypervolume
import matplotlib.pyplot as plt

### Parameters
assigning_problem = False
random_mode = 1 # 1 - random data, 2 - epsilon MOEA

num_runs = 10

save_path = "C:\\Users\\rosha\\Documents\\SEAK Lab Github\\VASSAR\\VASSAR_exec_heur\\results\\Operator Metrics\\" # for laptop
#save_path = "C:\\SEAK Lab Github\\VASSAR\\VASSAR_exec_heur\\results" # for workstation

filepath_prob = "Assigning\\"
if (not assigning_problem):
    filepath_prob = "Partitioning\\"
    
filepath_mode = "Random\\"
if assigning_problem:
    filename_mode = "random_assign_operator_index_"
else:
    filename_mode = "random_partition_operator_index_"
#filename_finalpop = ""    

if (random_mode == 2):
    filepath_mode = "Epsilon MOEA\\"
    if assigning_problem:
        filename_mode = "EpsilonMOEA_emoea_ClimateCentric_assign_operator_data_"
    else:
        filename_mode = "EpsilonMOEA_emoea_ClimateCentric_partition_operator_data_"
    filename_moea = filename_mode
    
    #if not final_pop_only:
        #filename_finalpop = "_fullpop"
        
filepath_full = save_path + filepath_prob + filepath_mode 

### Useful functions
def find_closest_index(val, search_list):
    val_diff = np.array(search_list) - val
    closest_index = np.argmin(np.abs(val_diff))
    return closest_index

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

def read_csv_run(filepath, filename, run_mode, nfe_thresh):
    full_filename = filepath + filename 
    
    with open(full_filename,newline='') as csvfile:
        data = [row for row in csv.reader(csvfile)]
        
        if (run_mode == 2):
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
        if (run_mode == 2):
            buffer = 1
        
        valid_count = 0
        for i in range(len(data)-1):
            
            if (run_mode == 2):
                nfe = list(map(float,data[i+1][buffer]))
            
            data_float = list(map(float,data[i+1][buffer+1:buffer+2]))
            data_float_instrdc = list(map(float,data[i+1][buffer+4:buffer+5]))
            data_float_instrorb = list(map(float,data[i+1][buffer+7:buffer+8]))
            data_float_interinstr = list(map(float,data[i+1][buffer+9:buffer+10]))
            data_float_packeff = list(map(float,data[i+1][buffer+11:buffer+12]))
            data_float_spmass = list(map(float,data[i+1][buffer+13:buffer+14]))
            data_float_instrsyn = list(map(float,data[i+1][buffer+15:buffer+16]))
            
            if (run_mode == 2):
                data_all = nfes + data_float + data_float_instrdc + data_float_instrorb + data_float_interinstr + data_float_packeff + data_float_spmass + data_float_instrsyn
            else:
                data_all = data_float + data_float_instrdc + data_float_instrorb + data_float_interinstr + data_float_packeff + data_float_spmass + data_float_instrsyn
            
            if (any(np.isnan(np.array(data_all))) or any(np.isinf(np.array(data_all)))):
                continue
            
            if (run_mode == 2):
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
            
    if (run_mode == 2):
        sort_indices = np.argsort(nfes)
        nfes_sorted = list(nfes[sort_indices])
        
        science_sorted = list(science[sort_indices])
        cost_sorted = list(cost[sort_indices])
        
        science_instrdc_sorted = list(science_instrdc[sort_indices])
        cost_instrdc_sorted = list(cost_instrdc[sort_indices])
        
        science_instrorb_sorted = list(science_instrorb[sort_indices])
        cost_instrorb_sorted = list(cost_instrorb[sort_indices])
        
        science_interinstr_sorted = list(science_interinstr[sort_indices])
        cost_interinstr_sorted = list(cost_interinstr[sort_indices])
        
        science_packeff_sorted = list(science_packeff[sort_indices])
        cost_packeff_sorted = list(cost_packeff[sort_indices])
        
        science_spmass_sorted = list(science_spmass[sort_indices])
        cost_spmass_sorted = list(cost_spmass[sort_indices])
        
        science_instrsyn_sorted = list(science_instrsyn[sort_indices])
        cost_instrsyn_sorted = list(cost_instrsyn[sort_indices])
        
        nfe_index = find_closest_index(nfe_thresh, nfes_sorted)
        
        objs = np.column_stack((science_sorted[:nfe_index], cost_sorted[:nfe_index]))
        objs_instrdc = np.column_stack((science_instrdc_sorted[:nfe_index], cost_instrdc_sorted[:nfe_index]))
        objs_instrorb = np.column_stack((science_instrorb_sorted[:nfe_index], cost_instrorb_sorted[:nfe_index]))
        objs_interinstr = np.column_stack((science_interinstr_sorted[:nfe_index], cost_interinstr_sorted[:nfe_index]))
        objs_packeff = np.column_stack((science_packeff_sorted[:nfe_index], cost_packeff_sorted[:nfe_index]))
        objs_spmass = np.column_stack((science_spmass_sorted[:nfe_index], cost_spmass_sorted[:nfe_index]))
        objs_instrsyn = np.column_stack((science_instrsyn_sorted[:nfe_index], cost_instrsyn_sorted[:nfe_index]))
    
    else:
        
        objs = np.column_stack((science, cost))
        objs_instrdc = np.column_stack((science_instrdc, cost_instrdc))
        objs_instrorb = np.column_stack((science_instrorb, cost_instrorb))
        objs_interinstr = np.column_stack((science_interinstr, cost_interinstr))
        objs_packeff = np.column_stack((science_packeff, cost_packeff))
        objs_spmass = np.column_stack((science_spmass, cost_spmass))
        objs_instrsyn = np.column_stack((science_instrsyn, cost_instrsyn))
    
    data = {}
    data['original'] = objs
    data['instrdc'] = objs_instrdc
    data['instrorb'] = objs_instrorb
    data['interinstr'] = objs_interinstr
    data['packeff'] = objs_packeff
    data['spmass'] = objs_spmass
    data['instrsyn'] = objs_instrsyn
    
    return data

### Operation
nfe_max = 5000 # for eps MOEA
n_datapoints = 10

nfes_array = np.linspace(0, nfe_max, n_datapoints+1)