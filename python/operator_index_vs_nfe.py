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
#random_mode = 1 # 1 - random data, 2 - epsilon MOEA

num_runs = 10

### Useful functions
def get_fileloc(assign_prob, rand_mode, run_num):
    #save_path = "C:\\Users\\rosha\\Documents\\SEAK Lab Github\\VASSAR\\VASSAR_exec_heur\\results\\Operator Metrics\\" # for laptop
    save_path = "C:\\SEAK Lab\\SEAK Lab Github\\VASSAR\\VASSAR_exec_heur\\results\\Operator Metrics\\" # for workstation

    filepath_prob = "Assigning\\"
    if (not assign_prob):
        filepath_prob = "Partitioning\\"

    if rand_mode == 1: 
        filepath_mode = 'Random\\'
        if assign_prob:
            filename_init = "random_assign_operator_index_"
        else:
            filename_init = "random_partition_operator_index_"   
        filename_mode = ''
    elif rand_mode == 2:
        if assigning_problem:
            filename_init = "EpsilonMOEA_emoea_ClimateCentric_assign_operator_data_"
        else:
            filename_init = "EpsilonMOEA_emoea_ClimateCentric_partition_operator_data_"
        filepath_mode = 'Epsilon MOEA\\'
        filename_mode = '_fullpop'

    filepath_full = save_path + filepath_prob + filepath_mode + filename_init + str(run_num) + filename_mode + '.csv'
    
    return filepath_full

def find_closest_index(val, search_list):
    val_diff = np.array(search_list) - val
    closest_index = np.argmin(np.abs(val_diff))
    return closest_index

def find_last_index(val,search_list):
    if val in search_list:
        idx = len(search_list) - search_list[::-1].index(val) - 1
    else:
        idx = 0
    return idx

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

def read_csv_run(full_filename, run_mode, nfe_thresh):
    
    with open(full_filename, newline='') as csvfile:
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
                nfe = int(data[i+1][buffer-1])
            
            data_float = list(map(float,data[i+1][buffer+1:buffer+3]))
            data_float_instrdc = list(map(float,data[i+1][buffer+4:buffer+6]))
            data_float_instrorb = list(map(float,data[i+1][buffer+7:buffer+9]))
            data_float_interinstr = list(map(float,data[i+1][buffer+10:buffer+12]))
            data_float_packeff = list(map(float,data[i+1][buffer+13:buffer+15]))
            data_float_spmass = list(map(float,data[i+1][buffer+16:buffer+18]))
            data_float_instrsyn = list(map(float,data[i+1][buffer+19:buffer+21]))
            
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
            
            valid_count += 1
            
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
        
        nfe_index = find_last_index(nfe_thresh, nfes_sorted)
        
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

def get_combined_data(data_rand, data_moea, assign_prob, rand_mode):
    
    if rand_mode == 1:
        objs_combined = data_rand['original'].tolist()
        objs_combined_instrdc = data_rand['instrdc'].tolist()
        objs_combined_instrorb = data_rand['instrorb'].tolist()
        objs_combined_interinstr = data_rand['interinstr'].tolist()
        objs_combined_packeff = data_rand['packeff'].tolist()
        objs_combined_spmass = data_rand['spmass'].tolist()
        objs_combined_instrsyn = data_rand['instrsyn'].tolist()
        
    else:
        objs_combined = data_rand['original'].tolist() + data_moea['original'].tolist()
        objs_combined_instrdc = data_rand['instrdc'].tolist() + data_moea['instrdc'].tolist()
        objs_combined_instrorb = data_rand['instrorb'].tolist() + data_moea['instrorb'].tolist()
        objs_combined_interinstr = data_rand['interinstr'].tolist() + data_moea['interinstr'].tolist()
        objs_combined_packeff = data_rand['packeff'].tolist() + data_moea['packeff'].tolist()
        objs_combined_spmass = data_rand['spmass'].tolist() + data_moea['spmass'].tolist()
        objs_combined_instrsyn = data_rand['instrsyn'].tolist() + data_moea['instrsyn'].tolist()
        
    # Compute Pareto Fronts
    pf_combined = compute_pareto_front(objs_combined)
    pf_combined_instrdc = compute_pareto_front(objs_combined_instrdc)
    pf_combined_instrorb = compute_pareto_front(objs_combined_instrorb)
    pf_combined_interinstr = compute_pareto_front(objs_combined_interinstr)
    pf_combined_packeff = compute_pareto_front(objs_combined_packeff)
    pf_combined_spmass = compute_pareto_front(objs_combined_spmass)
    pf_combined_instrsyn = compute_pareto_front(objs_combined_instrsyn)
    
    data_comb = {}
    data_comb['n_designs'] = len(objs_combined)
    
    data_comb['pf_original'] = pf_combined
    data_comb['pf_instrdc'] = pf_combined_instrdc
    data_comb['pf_instrorb'] = pf_combined_instrorb
    data_comb['pf_interinstr'] = pf_combined_interinstr
    data_comb['pf_packeff'] = pf_combined_packeff
    data_comb['pf_spmass'] = pf_combined_spmass
    data_comb['pf_instrsyn'] = pf_combined_instrsyn
    
    return data_comb

def get_obj_bounds_run(data_comb):
    
    pf = data_comb['pf_original']
    pf_instrdc = data_comb['pf_instrdc']
    pf_instrorb = data_comb['pf_instrorb']
    pf_interinstr = data_comb['pf_interinstr']
    pf_packeff = data_comb['pf_packeff']
    pf_spmass = data_comb['pf_spmass']
    pf_instrsyn = data_comb['pf_instrsyn']       
    
    pf_objs1 = [x[0] for x in pf]
    pf_objs2 = [x[1] for x in pf]
    
    pf_instrdc_objs1 = [x[0] for x in pf_instrdc]
    pf_instrdc_objs2 = [x[1] for x in pf_instrdc]
    
    pf_instrorb_objs1 = [x[0] for x in pf_instrorb]
    pf_instrorb_objs2 = [x[1] for x in pf_instrorb]
    
    pf_interinstr_objs1 = [x[0] for x in pf_interinstr]
    pf_interinstr_objs2 = [x[1] for x in pf_interinstr]
    
    pf_packeff_objs1 = [x[0] for x in pf_packeff]
    pf_packeff_objs2 = [x[1] for x in pf_packeff]
    
    pf_spmass_objs1 = [x[0] for x in pf_spmass]
    pf_spmass_objs2 = [x[1] for x in pf_spmass]
    
    pf_instrsyn_objs1 = [x[0] for x in pf_instrsyn]
    pf_instrsyn_objs2 = [x[1] for x in pf_instrsyn]
    
    objs1_max = np.amax([np.amax(pf_objs1), np.amax(pf_instrdc_objs1), np.amax(pf_instrorb_objs1), np.amax(pf_interinstr_objs1), np.amax(pf_packeff_objs1), np.amax(pf_spmass_objs1), np.amax(pf_instrsyn_objs1)])
    objs1_min = np.amin([np.amin(pf_objs1), np.amin(pf_instrdc_objs1), np.amin(pf_instrorb_objs1), np.amin(pf_interinstr_objs1), np.amin(pf_packeff_objs1), np.amin(pf_spmass_objs1), np.amin(pf_instrsyn_objs1)])
    objs2_max = np.amax([np.amax(pf_objs2), np.amax(pf_instrdc_objs2), np.amax(pf_instrorb_objs2), np.amax(pf_interinstr_objs2), np.amax(pf_packeff_objs2), np.amax(pf_spmass_objs2), np.amax(pf_instrsyn_objs2)])
    objs2_min = np.amin([np.amin(pf_objs2), np.amin(pf_instrdc_objs2), np.amin(pf_instrorb_objs2), np.amin(pf_interinstr_objs2), np.amin(pf_packeff_objs2), np.amin(pf_spmass_objs2), np.amin(pf_instrsyn_objs2)])
    
    return [objs1_max, objs1_min, objs2_max, objs2_min]

def compute_I_hv(pf_h, pf):
    return compute_hv(pf_h) - compute_hv(pf)
    
def compute_indices(data_comb, assign_prob, obj_bounds_all): # to be called after combining data
    pf = data_comb['pf_original']
    pf_instrdc = data_comb['pf_instrdc']
    pf_instrorb = data_comb['pf_instrorb']
    pf_interinstr = data_comb['pf_interinstr']
    pf_packeff = data_comb['pf_packeff']
    pf_spmass = data_comb['pf_spmass']
    pf_instrsyn = data_comb['pf_instrsyn']
    
    #n_total = data_comb['n_designs']
    
    # Normalize pareto fronts
    pf_objs1_norm = [(x[0] - obj_bounds_all[1])/(obj_bounds_all[0] - obj_bounds_all[1]) for x in pf]
    pf_objs1_instrdc_norm = [(x[0] - obj_bounds_all[1])/(obj_bounds_all[0] - obj_bounds_all[1]) for x in pf_instrdc]
    pf_objs1_instrorb_norm = [(x[0] - obj_bounds_all[1])/(obj_bounds_all[0] - obj_bounds_all[1]) for x in pf_instrorb]
    pf_objs1_interinstr_norm = [(x[0] - obj_bounds_all[1])/(obj_bounds_all[0] - obj_bounds_all[1]) for x in pf_interinstr]
    pf_objs1_packeff_norm = [(x[0] - obj_bounds_all[1])/(obj_bounds_all[0] - obj_bounds_all[1]) for x in pf_packeff]
    pf_objs1_spmass_norm = [(x[0] - obj_bounds_all[1])/(obj_bounds_all[0] - obj_bounds_all[1]) for x in pf_spmass]
    pf_objs1_instrsyn_norm = [(x[0] - obj_bounds_all[1])/(obj_bounds_all[0] - obj_bounds_all[1]) for x in pf_instrsyn]
    
    pf_objs2_norm = [(x[1] - obj_bounds_all[3])/(obj_bounds_all[2] - obj_bounds_all[3]) for x in pf]
    pf_objs2_instrdc_norm = [(x[1] - obj_bounds_all[3])/(obj_bounds_all[2] - obj_bounds_all[3]) for x in pf_instrdc]
    pf_objs2_instrorb_norm = [(x[1] - obj_bounds_all[3])/(obj_bounds_all[2] - obj_bounds_all[3]) for x in pf_instrorb]
    pf_objs2_interinstr_norm = [(x[1] - obj_bounds_all[3])/(obj_bounds_all[2] - obj_bounds_all[3]) for x in pf_interinstr]
    pf_objs2_packeff_norm = [(x[1] - obj_bounds_all[3])/(obj_bounds_all[2] - obj_bounds_all[3]) for x in pf_packeff]
    pf_objs2_spmass_norm = [(x[1] - obj_bounds_all[3])/(obj_bounds_all[2] - obj_bounds_all[3]) for x in pf_spmass]
    pf_objs2_instrsyn_norm = [(x[1] - obj_bounds_all[3])/(obj_bounds_all[2] - obj_bounds_all[3]) for x in pf_instrsyn]
    
    pf_norm = np.column_stack((pf_objs1_norm, pf_objs2_norm))
    pf_instrdc_norm = np.column_stack((pf_objs1_instrdc_norm, pf_objs2_instrdc_norm))
    pf_instrorb_norm = np.column_stack((pf_objs1_instrorb_norm, pf_objs2_instrorb_norm))
    pf_interinstr_norm = np.column_stack((pf_objs1_interinstr_norm, pf_objs2_interinstr_norm))
    pf_packeff_norm = np.column_stack((pf_objs1_packeff_norm, pf_objs2_packeff_norm))
    pf_spmass_norm = np.column_stack((pf_objs1_spmass_norm, pf_objs2_spmass_norm))
    pf_instrsyn_norm = np.column_stack((pf_objs1_instrsyn_norm, pf_objs2_instrsyn_norm))
    
    # Compute indices based on objective minimization
    I_instrdc = compute_I_hv(pf_instrdc_norm, pf_norm) 
    I_instrorb = compute_I_hv(pf_instrorb_norm, pf_norm) 
    I_interinstr = compute_I_hv(pf_interinstr_norm, pf_norm)
    I_packeff = compute_I_hv(pf_packeff_norm, pf_norm) 
    I_spmass = compute_I_hv(pf_spmass_norm, pf_norm) 
    I_instrsyn = compute_I_hv(pf_instrsyn_norm, pf_norm)
    
    return I_instrdc, I_instrorb, I_interinstr, I_packeff, I_spmass, I_instrsyn

def compute_indices_allruns(n_runs, assign_prob, nfe_threshold):
    I_instrdc_allruns = np.zeros(n_runs)
    I_instrorb_allruns = np.zeros(n_runs)
    I_interinstr_allruns = np.zeros(n_runs)
    I_packeff_allruns = np.zeros(n_runs)
    I_spmass_allruns = np.zeros(n_runs)
    I_instrsyn_allruns = np.zeros(n_runs)
    
    ### Compute Pareto Fronts for all runs first
    objs1_max_allruns = np.zeros(n_runs)
    objs1_min_allruns = np.zeros(n_runs)
    objs2_max_allruns = np.zeros(n_runs)
    objs2_min_allruns = np.zeros(n_runs)
    data_comb_allruns = {}
    for i in range(n_runs):
        filepath_rand = get_fileloc(assigning_problem, 1, i)
        data_rand_i = read_csv_run(filepath_rand, 1, nfe_threshold)
        
        filepath_moea = get_fileloc(assigning_problem, 2, i)
        data_moea_i = read_csv_run(filepath_moea, 2, nfe_threshold)

        data_comb_i = get_combined_data(data_rand_i, data_moea_i, assigning_problem, 2)
        
        data_comb_allruns['run'+str(i)] = data_comb_i
        
        obj_bounds_i = get_obj_bounds_run(data_comb_i)
        objs1_max_allruns[i] = obj_bounds_i[0]
        objs1_min_allruns[i] = obj_bounds_i[1]
        objs2_max_allruns[i] = obj_bounds_i[2]
        objs2_min_allruns[i] = obj_bounds_i[3]
    
    ### Compute overall objective bounds 
    objs1_max_overall = np.amax(objs1_max_allruns)
    objs1_min_overall = np.amin(objs1_min_allruns)
    objs2_max_overall = np.amax(objs2_max_allruns)
    objs2_min_overall = np.amin(objs2_min_allruns)
    obj_bounds_overall = [objs1_max_overall, objs1_min_overall, objs2_max_overall, objs2_min_overall] 

    ### Compute indices
    for i in range(n_runs):
        data_comb_i = data_comb_allruns['run'+str(i)]
        I_instrdci, I_instrorbi, I_interinstri, I_packeffi, I_spmassi, I_instrsyni = compute_indices(data_comb_i, assigning_problem, obj_bounds_overall)
        I_instrdc_allruns[i] = I_instrdci
        I_instrorb_allruns[i] = I_instrorbi
        I_interinstr_allruns[i] = I_interinstri
        I_packeff_allruns[i] = I_packeffi
        I_spmass_allruns[i] = I_spmassi
        I_instrsyn_allruns[i] = I_instrsyni
    
    return I_instrdc_allruns, I_instrorb_allruns, I_interinstr_allruns, I_packeff_allruns, I_spmass_allruns, I_instrsyn_allruns    

### Operation
nfe_max = 5000 # for eps MOEA
n_datapoints = 10

nfes_array = np.linspace(0, nfe_max, n_datapoints+1)

nfes_array_f = np.linspace(0, nfe_max, n_datapoints+1)
nfes_array = [int(x) for x in nfes_array_f]

I_instrdc_allnfes = []
I_instrorb_allnfes = []
I_interinstr_allnfes = []
I_packeff_allnfes = []
I_spmass_allnfes = []
I_instrsyn_allnfes = []

for nfe in nfes_array:
    I_instrdc_nfe, I_instrorb_nfe, I_interinstr_nfe, I_packeff_nfe, I_spmass_nfe, I_instrsyn_nfe = compute_indices_allruns(num_runs, assigning_problem, nfe)
    I_instrdc_allnfes.append(I_instrdc_nfe)
    I_instrorb_allnfes.append(I_instrorb_nfe)
    I_interinstr_allnfes.append(I_interinstr_nfe)
    I_packeff_allnfes.append(I_packeff_nfe)
    I_spmass_allnfes.append(I_spmass_nfe)
    I_instrsyn_allnfes.append(I_instrsyn_nfe)
    
### Plotting
nfe_idx_ticks = np.linspace(1, len(nfes_array), len(nfes_array))

plt.figure()
plt.boxplot(I_instrdc_allnfes)
plt.xticks(nfe_idx_ticks, nfes_array, fontsize=8)
plt.xlabel('NFE')
plt.ylabel('Heuristic Index', fontsize=8)
plt.title('Duty Cycle Viol.', fontsize=8)

plt.figure()
plt.boxplot(I_instrorb_allnfes)
plt.xticks(nfe_idx_ticks, nfes_array, fontsize=8)
plt.xlabel('NFE')
plt.ylabel('Heuristic Index', fontsize=8)
plt.title('Instr. Orbit Rel. Viol.', fontsize=8)

plt.figure()
plt.boxplot(I_interinstr_allnfes)
plt.xticks(nfe_idx_ticks, nfes_array, fontsize=8)
plt.xlabel('NFE')
plt.ylabel('Heuristic Index', fontsize=8)
plt.title('Interference Viol.', fontsize=8)

plt.figure()
plt.boxplot(I_packeff_allnfes)
plt.xticks(nfe_idx_ticks, nfes_array, fontsize=8)
plt.xlabel('NFE')
plt.ylabel('Heuristic Index', fontsize=8)
plt.title('Pack. Eff. Viol.', fontsize=8)

plt.figure()
plt.boxplot(I_spmass_allnfes)
plt.xticks(nfe_idx_ticks, nfes_array, fontsize=8)
plt.xlabel('NFE')
plt.ylabel('Heuristic Index', fontsize=8)
plt.title('Spacecraft Mass Viol.', fontsize=8)

plt.figure()
plt.boxplot(I_instrsyn_allnfes)
plt.xticks(nfe_idx_ticks, nfes_array, fontsize=8)
plt.xlabel('NFE')
plt.ylabel('Heuristic Index', fontsize=8)
plt.title('Synergy Viol.', fontsize=8)
