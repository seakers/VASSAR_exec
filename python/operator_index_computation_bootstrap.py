# -*- coding: utf-8 -*-
"""
Operator Index for satellite design problems with Bootstrapping

@author: roshan94
"""
import csv
import numpy as np
from pygmo import hypervolume

### Parameters
assigning_problem = True
random_mode = 1 # 1 - only random data, 2 - only epsilon MOEA, 3 - both (always 1)
#final_pop_only = True # only specific to epsilon MOEA

#num_runs = 10

#num_archs = 300 # only specific to random data 

#save_path = "C:\\Users\\rosha\\Documents\\SEAK Lab Github\\VASSAR\\VASSAR_exec_heur\\results\\Operator Metrics\\" # for laptop
save_path = "C:\\SEAK Lab\\SEAK Lab Github\\VASSAR\\VASSAR_exec_heur\\results\\Operator Metrics\\" # for workstation

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
    
    #if not final_pop_only:
        #filename_finalpop = "_fullpop"
        
filepath_full = save_path + filepath_prob + filepath_mode 

obj_weights = [1, 1]
if assigning_problem:
    obj_weights = [0.425, 2.5e4] # normalization weights for objectives
    
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

def read_csv_run(filepath, filename, run_mode, assign_prob):
    full_filename = filepath + filename 
    
    with open(full_filename,newline='') as csvfile:
        data = [row for row in csv.reader(csvfile)]

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
        
        if assign_prob:
            science_instrcount = np.zeros(len(data)-1)
            cost_instrcount = np.zeros(len(data)-1)
        
        valid_count = 0
        for i in range(len(data)-1):
            
            data_float = list(map(float,data[i+1][1:3]))
            data_float_instrdc = list(map(float,data[i+1][4:6]))
            data_float_instrorb = list(map(float,data[i+1][7:9]))
            data_float_interinstr = list(map(float,data[i+1][10:12]))
            data_float_packeff = list(map(float,data[i+1][13:15]))
            data_float_spmass = list(map(float,data[i+1][16:18]))
            data_float_instrsyn = list(map(float,data[i+1][19:21]))
            if assign_prob:
                data_float_instrcount = list(map(float,data[i+1][22:24]))
            
            data_all = data_float + data_float_instrdc + data_float_instrorb + data_float_interinstr + data_float_packeff + data_float_spmass + data_float_instrsyn
            if assign_prob:
                data_all = data_all + data_float_instrcount
            
            if (any(np.isnan(np.array(data_all))) or any(np.isinf(np.array(data_all)))):
                continue
            
            science[valid_count] = data_float[0]*obj_weights[0]
            cost[valid_count] = data_float[1]*obj_weights[1]
            
            science_instrdc[valid_count] = data_float_instrdc[0]*obj_weights[0]
            cost_instrdc[valid_count] = data_float_instrdc[1]*obj_weights[1]
            
            science_instrorb[valid_count] = data_float_instrorb[0]*obj_weights[0]
            cost_instrorb[valid_count] = data_float_instrorb[1]*obj_weights[1]
            
            science_interinstr[valid_count] = data_float_interinstr[0]*obj_weights[0]
            cost_interinstr[valid_count] = data_float_interinstr[1]*obj_weights[1]
            
            science_packeff[valid_count] = data_float_packeff[0]*obj_weights[0]
            cost_packeff[valid_count] = data_float_packeff[1]*obj_weights[1]
            
            science_spmass[valid_count] = data_float_spmass[0]*obj_weights[0]
            cost_spmass[valid_count] = data_float_spmass[1]*obj_weights[1]
            
            science_instrsyn[valid_count] = data_float_instrsyn[0]*obj_weights[0]
            cost_instrsyn[valid_count] = data_float_instrsyn[1]*obj_weights[1]
            
            if assign_prob:
                science_instrcount[valid_count] = data_float_instrcount[0]*obj_weights[0]
                cost_instrcount[valid_count] = data_float_instrcount[1]*obj_weights[1]
            
            valid_count += 1
            
    objs = np.column_stack((science, cost))
    objs_instrdc = np.column_stack((science_instrdc, cost_instrdc))
    objs_instrorb = np.column_stack((science_instrorb, cost_instrorb))
    objs_interinstr = np.column_stack((science_interinstr, cost_interinstr))
    objs_packeff = np.column_stack((science_packeff, cost_packeff))
    objs_spmass = np.column_stack((science_spmass, cost_spmass))
    objs_instrsyn = np.column_stack((science_instrsyn, cost_instrsyn))
    if assign_prob:
        objs_instrcount = np.column_stack((science_instrcount, cost_instrcount))
        
    data = {}
    data['original'] = objs
    data['instrdc'] = objs_instrdc
    data['instrorb'] = objs_instrorb
    data['interinstr'] = objs_interinstr
    data['packeff'] = objs_packeff
    data['spmass'] = objs_spmass
    data['instrsyn'] = objs_instrsyn
    if assign_prob:
        data['instrcount'] = objs_instrcount
    
    return data

def get_data_subset(data, dataset_inds, assign_prob):
    dataset = {}
    
    dataset['original'] = data['original'][dataset_inds,:]
    dataset['instrdc'] = data['instrdc'][dataset_inds,:]
    dataset['instrorb'] = data['instrorb'][dataset_inds,:]
    dataset['interinstr'] = data['interinstr'][dataset_inds,:]
    dataset['packeff'] = data['packeff'][dataset_inds,:]
    dataset['spmass'] = data['spmass'][dataset_inds,:]
    dataset['instrsyn'] = data['instrsyn'][dataset_inds,:]
    if assign_prob:
        dataset['instrcount'] = data['instrcount'][dataset_inds,:]
        
    return dataset

def compute_indices(data, assign_prob):
    
    data_orig = data['original']
    data_instrdc = data['instrdc']
    data_instrorb = data['instrorb']
    data_interinstr = data['interinstr']
    data_packeff = data['packeff']
    data_spmass = data['spmass']
    data_instrsyn = data['instrsyn']
    if assigning_problem:
        data_instrcount = data['instrcount']
        
    ## Compute Pareto Fronts
    pf_orig = compute_pareto_front(data_orig)
    pf_instrdc = compute_pareto_front(data_instrdc)
    pf_instrorb = compute_pareto_front(data_instrorb)
    pf_interinstr = compute_pareto_front(data_interinstr)
    pf_packeff = compute_pareto_front(data_packeff)
    pf_spmass = compute_pareto_front(data_spmass)
    pf_instrsyn = compute_pareto_front(data_instrsyn)
    if assigning_problem:
        pf_instrcount = compute_pareto_front(data_instrcount)
        
    ## Obtain bounds for normalization
    science_pf_orig = [x[0] for x in pf_orig]
    cost_pf_orig = [x[1] for x in pf_orig]
    
    science_pf_instrdc = [x[0] for x in pf_instrdc]
    cost_pf_instrdc = [x[1] for x in pf_instrdc]
    
    science_pf_instrorb = [x[0] for x in pf_instrorb]
    cost_pf_instrorb = [x[1] for x in pf_instrorb]
    
    science_pf_interinstr = [x[0] for x in pf_interinstr]
    cost_pf_interinstr = [x[1] for x in pf_interinstr]
    
    science_pf_packeff = [x[0] for x in pf_packeff]
    cost_pf_packeff = [x[1] for x in pf_packeff]
    
    science_pf_spmass = [x[0] for x in pf_spmass]
    cost_pf_spmass = [x[1] for x in pf_spmass]
    
    science_pf_instrsyn = [x[0] for x in pf_instrsyn]
    cost_pf_instrsyn = [x[1] for x in pf_instrsyn]
    
    if assigning_problem:
        science_pf_instrcount = [x[0] for x in pf_instrcount]
        cost_pf_instrcount = [x[1] for x in pf_instrcount]
        
        science_norm_max = np.amax([np.amax(science_pf_orig), np.amax(science_pf_instrdc), np.amax(science_pf_instrorb), np.amax(science_pf_interinstr), np.amax(science_pf_packeff), np.amax(science_pf_spmass), np.amax(science_pf_instrsyn), np.amax(science_pf_instrcount)])
        cost_norm_max = np.amax([np.amax(cost_pf_orig), np.amax(cost_pf_instrdc), np.amax(cost_pf_instrorb), np.amax(cost_pf_interinstr), np.amax(cost_pf_packeff), np.amax(cost_pf_spmass), np.amax(cost_pf_instrsyn), np.amax(cost_pf_instrcount)])
            
        science_norm_min = np.amin([np.amin(science_pf_orig), np.amin(science_pf_instrdc), np.amin(science_pf_instrorb), np.amin(science_pf_interinstr), np.amin(science_pf_packeff), np.amin(science_pf_spmass), np.amin(science_pf_instrsyn), np.amin(science_pf_instrcount)])
        cost_norm_min = np.amin([np.amin(cost_pf_orig), np.amin(cost_pf_instrdc), np.amin(cost_pf_instrorb), np.amin(cost_pf_interinstr), np.amin(cost_pf_packeff), np.amin(cost_pf_spmass), np.amin(cost_pf_instrsyn), np.amin(cost_pf_instrcount)])
        
    else:
        science_norm_max = np.amax([np.amax(science_pf_orig), np.amax(science_pf_instrdc), np.amax(science_pf_instrorb), np.amax(science_pf_interinstr), np.amax(science_pf_packeff), np.amax(science_pf_spmass), np.amax(science_pf_instrsyn)])
        cost_norm_max = np.amax([np.amax(cost_pf_orig), np.amax(cost_pf_instrdc), np.amax(cost_pf_instrorb), np.amax(cost_pf_interinstr), np.amax(cost_pf_packeff), np.amax(cost_pf_spmass), np.amax(cost_pf_instrsyn)])
            
        science_norm_min = np.amin([np.amin(science_pf_orig), np.amin(science_pf_instrdc), np.amin(science_pf_instrorb), np.amin(science_pf_interinstr), np.amin(science_pf_packeff), np.amin(science_pf_spmass), np.amin(science_pf_instrsyn)])
        cost_norm_min = np.amin([np.amin(cost_pf_orig), np.amin(cost_pf_instrdc), np.amin(cost_pf_instrorb), np.amin(cost_pf_interinstr), np.amin(cost_pf_packeff), np.amin(cost_pf_spmass), np.amin(cost_pf_instrsyn)])
        
    ## Obtain normalized Pareto Fronts
    science_pf_norm_orig = [((x[0] - science_norm_min)/(science_norm_max - science_norm_min)) for x in pf_orig]
    cost_pf_norm_orig = [((x[1] - cost_norm_min)/(cost_norm_max - cost_norm_min)) for x in pf_orig]
    
    science_pf_norm_instrdc = [((x[0] - science_norm_min)/(science_norm_max - science_norm_min)) for x in pf_instrdc]
    cost_pf_norm_instrdc = [((x[1] - cost_norm_min)/(cost_norm_max - cost_norm_min)) for x in pf_instrdc]
    
    science_pf_norm_instrorb = [((x[0] - science_norm_min)/(science_norm_max - science_norm_min)) for x in pf_instrorb]
    cost_pf_norm_instrorb = [((x[1] - cost_norm_min)/(cost_norm_max - cost_norm_min)) for x in pf_instrorb]
    
    science_pf_norm_interinstr = [((x[0] - science_norm_min)/(science_norm_max - science_norm_min)) for x in pf_interinstr]
    cost_pf_norm_interinstr = [((x[1] - cost_norm_min)/(cost_norm_max - cost_norm_min)) for x in pf_interinstr]
    
    science_pf_norm_packeff = [((x[0] - science_norm_min)/(science_norm_max - science_norm_min)) for x in pf_packeff]
    cost_pf_norm_packeff = [((x[1] - cost_norm_min)/(cost_norm_max - cost_norm_min)) for x in pf_packeff]
    
    science_pf_norm_spmass = [((x[0] - science_norm_min)/(science_norm_max - science_norm_min)) for x in pf_spmass]
    cost_pf_norm_spmass = [((x[1] - cost_norm_min)/(cost_norm_max - cost_norm_min)) for x in pf_spmass]
    
    science_pf_norm_instrsyn = [((x[0] - science_norm_min)/(science_norm_max - science_norm_min)) for x in pf_instrsyn]
    cost_pf_norm_instrsyn = [((x[1] - cost_norm_min)/(cost_norm_max - cost_norm_min)) for x in pf_instrsyn]
    
    if assigning_problem:
        science_pf_norm_instrcount = [((x[0] - science_norm_min)/(science_norm_max - science_norm_min)) for x in pf_instrcount]
        cost_pf_norm_instrcount = [((x[1] - cost_norm_min)/(cost_norm_max - cost_norm_min)) for x in pf_instrcount]
    
    pf_norm_orig = np.column_stack((science_pf_norm_orig, cost_pf_norm_orig))
    pf_norm_instrdc = np.column_stack((science_pf_norm_instrdc, cost_pf_norm_instrdc))
    pf_norm_instrorb = np.column_stack((science_pf_norm_instrorb, cost_pf_norm_instrorb))
    pf_norm_interinstr = np.column_stack((science_pf_norm_interinstr, cost_pf_norm_interinstr))
    pf_norm_packeff = np.column_stack((science_pf_norm_packeff, cost_pf_norm_packeff))
    pf_norm_spmass = np.column_stack((science_pf_norm_spmass, cost_pf_norm_spmass))
    pf_norm_instrsyn = np.column_stack((science_pf_norm_instrsyn, cost_pf_norm_instrsyn))
    if assigning_problem:
        pf_norm_instrcount = np.column_stack((science_pf_norm_instrcount, cost_pf_norm_instrcount))
    
    ## Compute Hypervolume and Indices
    hv_orig = compute_hv(pf_norm_orig)
    hv_instrdc = compute_hv(pf_norm_instrdc)
    hv_instrorb = compute_hv(pf_norm_instrorb)
    hv_interinstr = compute_hv(pf_norm_interinstr)
    hv_packeff = compute_hv(pf_norm_packeff)
    hv_spmass = compute_hv(pf_norm_spmass)
    hv_instrsyn = compute_hv(pf_norm_instrsyn)
    if assigning_problem:
        hv_instrcount = compute_hv(pf_norm_instrcount)
    
    I_instrdc = hv_instrdc - hv_orig
    I_instrorb = hv_instrorb - hv_orig
    I_interinstr = hv_interinstr - hv_orig
    I_packeff = hv_packeff - hv_orig
    I_spmass = hv_spmass - hv_orig
    I_instrsyn = hv_instrsyn - hv_orig
    if assigning_problem:
        I_instrcount = hv_instrcount - hv_orig
    else:
        I_instrcount = 0
        
    return [I_instrdc, I_instrorb, I_interinstr, I_packeff, I_spmass, I_instrsyn, I_instrcount]

#### OPERATION    
run_num = 0

filename_full = filename_mode + str(run_num) + '.csv'
data_all = read_csv_run(filepath_full, filename_full, random_mode, assigning_problem)

n_heurs = 6
if assigning_problem:
    n_heurs = 7

dataset_mode = 3 # 1 - bootstrapping, 2 - leave one out, 3 - equal division

if dataset_mode == 1:
    n_datasets = 30
    n_des = 200
elif dataset_mode == 2:
    n_datasets = 300
    n_des = 299
else:
    n_datasets = 10
    n_des = 30
    
### Compute Indices for each dataset
I_datasets = np.zeros((n_datasets, 7))
for i in range(n_datasets):
    
    if dataset_mode == 1:
        dataset_des_inds = np.random.randint(0, len(data_all['original']), n_des)
    elif dataset_mode == 2:
        dataset_des_inds_all = np.linspace(0, n_datasets-1, n_datasets)
        dataset_des_inds = np.delete(dataset_des_inds_all, i)
        dataset_des_inds = np.asarray(dataset_des_inds, dtype='int')
    else:
        dataset_des_inds = np.linspace(i*n_des, (i+1)*n_des-1, n_des)
        dataset_des_inds = np.asarray(dataset_des_inds, dtype='int')
        
    data_comb_run = get_data_subset(data_all, dataset_des_inds, assigning_problem)
    I_datasets[i,:] = compute_indices(data_comb_run, assigning_problem)
    
### Empirical probability of positive index
prob_heurs = np.zeros((n_heurs))
for i in range(n_heurs):
    prob_heurs[i] = len([x for x in I_datasets[:,i] if x > 0])/n_datasets