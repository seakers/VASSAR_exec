# -*- coding: utf-8 -*-
"""
Operator Index determination for satellite problems

@author: rosha
"""

import csv
import numpy as np
from pygmo import hypervolume

### Parameters
assigning_problem = True
random_mode = 1 # 1 - only random data, 2 - only epsilon MOEA, 3 - both
#final_pop_only = True # only specific to epsilon MOEA

num_runs = 10

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

def read_csv_run(filepath, filename, run_mode):
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
        
        valid_count = 0
        for i in range(len(data)-1):
            
            data_float = list(map(float,data[i+1][1:3]))
            data_float_instrdc = list(map(float,data[i+1][4:6]))
            data_float_instrorb = list(map(float,data[i+1][7:9]))
            data_float_interinstr = list(map(float,data[i+1][10:12]))
            data_float_packeff = list(map(float,data[i+1][13:15]))
            data_float_spmass = list(map(float,data[i+1][16:18]))
            data_float_instrsyn = list(map(float,data[i+1][19:21]))
            
            data_all = data_float + data_float_instrdc + data_float_instrorb + data_float_interinstr + data_float_packeff + data_float_spmass + data_float_instrsyn
            
            if (any(np.isnan(np.array(data_all))) or any(np.isinf(np.array(data_all)))):
                continue
            
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
I_instrdc = np.zeros(num_runs)
I_instrorb = np.zeros(num_runs)
I_interinstr = np.zeros(num_runs)
I_packeff = np.zeros(num_runs)
I_spmass = np.zeros(num_runs)
I_instrsyn = np.zeros(num_runs)

pfs_orig = {}
pfs_instrdc = {}
pfs_instrorb = {}
pfs_interinstr = {}
pfs_packeff = {}
pfs_spmass = {}
pfs_instrsyn = {}

science_max_allruns = np.zeros(num_runs)
cost_max_allruns = np.zeros(num_runs)
science_min_allruns = np.zeros(num_runs)
cost_min_allruns = np.zeros(num_runs)

for i in range(num_runs):
    filename_full = filename_mode + str(i) + '.csv'
    data_all = read_csv_run(filepath_full, filename_full, random_mode)
    
    if (random_mode == 3):
        if assigning_problem:
            filename_moea = "EpsilonMOEA_emoea_ClimateCentric_assign_operator_data_"
        else:
            filename_moea = "EpsilonMOEA_emoea_ClimateCentric_partition_operator_data_"
            
        filepath_moea_full = save_path + filepath_prob + "Epsilon MOEA\\"
            
        filename_moea_full = filename_moea + str(i) + '.csv'
        data_moea = read_csv_run(filepath_moea_full, filename_moea_full, random_mode)
    
        data_orig_combined = np.vstack((data_all['original'], data_moea['original']))
        data_instrdc_combined = np.vstack((data_all['instrdc'], data_moea['instrdc']))
        data_instrorb_combined = np.vstack((data_all['instrorb'], data_moea['instrorb']))
        data_interinstr_combined = np.vstack((data_all['interinstr'], data_moea['interinstr']))
        data_packeff_combined = np.vstack((data_all['packeff'], data_moea['packeff']))
        data_spmass_combined = np.vstack((data_all['spmass'], data_moea['spmass']))
        data_instrsyn_combined = np.vstack((data_all['instrsyn'], data_moea['instrsyn']))
        
    else:
        
        data_orig_combined = data_all['original']
        data_instrdc_combined = data_all['instrdc']
        data_instrorb_combined = data_all['instrorb']
        data_interinstr_combined = data_all['interinstr']
        data_packeff_combined = data_all['packeff']
        data_spmass_combined = data_all['spmass']
        data_instrsyn_combined = data_all['instrsyn']
    
    pf_orig = compute_pareto_front(data_orig_combined)
    pf_instrdc = compute_pareto_front(data_instrdc_combined)
    pf_instrorb = compute_pareto_front(data_instrorb_combined)
    pf_interinstr = compute_pareto_front(data_interinstr_combined)
    pf_packeff = compute_pareto_front(data_packeff_combined)
    pf_spmass = compute_pareto_front(data_spmass_combined)
    pf_instrsyn = compute_pareto_front(data_instrsyn_combined)
    
    pfs_orig['run'+str(i)] = pf_orig
    pfs_instrdc['run'+str(i)] = pf_instrdc
    pfs_instrorb['run'+str(i)] = pf_instrorb
    pfs_interinstr['run'+str(i)] = pf_interinstr
    pfs_packeff['run'+str(i)] = pf_packeff
    pfs_spmass['run'+str(i)] = pf_spmass
    pfs_instrsyn['run'+str(i)] = pf_instrsyn
    
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
    
    science_norm_max = np.amax([np.amax(science_pf_orig), np.amax(science_pf_instrdc), np.amax(science_pf_instrorb), np.amax(science_pf_interinstr), np.amax(science_pf_packeff), np.amax(science_pf_spmass), np.amax(science_pf_instrsyn)])
    cost_norm_max = np.amax([np.amax(cost_pf_orig), np.amax(cost_pf_instrdc), np.amax(cost_pf_instrorb), np.amax(cost_pf_interinstr), np.amax(cost_pf_packeff), np.amax(cost_pf_spmass), np.amax(cost_pf_instrsyn)])
        
    science_norm_min = np.amin([np.amin(science_pf_orig), np.amin(science_pf_instrdc), np.amin(science_pf_instrorb), np.amin(science_pf_interinstr), np.amin(science_pf_packeff), np.amin(science_pf_spmass), np.amin(science_pf_instrsyn)])
    cost_norm_min = np.amin([np.amin(cost_pf_orig), np.amin(cost_pf_instrdc), np.amin(cost_pf_instrorb), np.amin(cost_pf_interinstr), np.amin(cost_pf_packeff), np.amin(cost_pf_spmass), np.amin(cost_pf_instrsyn)])
            
    science_max_allruns[i] = science_norm_max
    cost_max_allruns[i] = cost_norm_max
    science_min_allruns[i] = science_norm_min
    cost_min_allruns[i] = cost_norm_min
    
science_norm_max_allruns = np.amax(science_max_allruns)
cost_norm_max_allruns = np.amax(cost_max_allruns)
science_norm_min_allruns = np.amin(science_min_allruns)
cost_norm_min_allruns = np.amin(cost_min_allruns)
            
for i in range(num_runs):
    pf_orig = pfs_orig['run'+str(i)]
    pf_instrdc = pfs_instrdc['run'+str(i)]
    pf_instrorb = pfs_instrorb['run'+str(i)]
    pf_interinstr = pfs_interinstr['run'+str(i)]
    pf_packeff = pfs_packeff['run'+str(i)]
    pf_spmass = pfs_spmass['run'+str(i)]
    pf_instrsyn = pfs_instrsyn['run'+str(i)]
    
    science_pf_norm_orig = [((x[0] - science_norm_min_allruns)/(science_norm_max_allruns - science_norm_min_allruns)) for x in pf_orig]
    cost_pf_norm_orig = [((x[1] - cost_norm_min_allruns)/(cost_norm_max_allruns - cost_norm_min_allruns)) for x in pf_orig]
    
    science_pf_norm_instrdc = [((x[0] - science_norm_min_allruns)/(science_norm_max_allruns - science_norm_min_allruns)) for x in pf_instrdc]
    cost_pf_norm_instrdc = [((x[1] - cost_norm_min_allruns)/(cost_norm_max_allruns - cost_norm_min_allruns)) for x in pf_instrdc]
    
    science_pf_norm_instrorb = [((x[0] - science_norm_min_allruns)/(science_norm_max_allruns - science_norm_min_allruns)) for x in pf_instrorb]
    cost_pf_norm_instrorb = [((x[1] - cost_norm_min_allruns)/(cost_norm_max_allruns - cost_norm_min_allruns)) for x in pf_instrorb]
    
    science_pf_norm_interinstr = [((x[0] - science_norm_min_allruns)/(science_norm_max_allruns - science_norm_min_allruns)) for x in pf_interinstr]
    cost_pf_norm_interinstr = [((x[1] - cost_norm_min_allruns)/(cost_norm_max_allruns - cost_norm_min_allruns)) for x in pf_interinstr]
    
    science_pf_norm_packeff = [((x[0] - science_norm_min_allruns)/(science_norm_max_allruns - science_norm_min_allruns)) for x in pf_packeff]
    cost_pf_norm_packeff = [((x[1] - cost_norm_min_allruns)/(cost_norm_max_allruns - cost_norm_min_allruns)) for x in pf_packeff]
    
    science_pf_norm_spmass = [((x[0] - science_norm_min_allruns)/(science_norm_max_allruns - science_norm_min_allruns)) for x in pf_spmass]
    cost_pf_norm_spmass = [((x[1] - cost_norm_min_allruns)/(cost_norm_max_allruns - cost_norm_min_allruns)) for x in pf_spmass]
    
    science_pf_norm_instrsyn = [((x[0] - science_norm_min_allruns)/(science_norm_max_allruns - science_norm_min_allruns)) for x in pf_instrsyn]
    cost_pf_norm_instrsyn = [((x[1] - cost_norm_min_allruns)/(cost_norm_max_allruns - cost_norm_min_allruns)) for x in pf_instrsyn]
    
    pf_norm_orig = np.column_stack((science_pf_norm_orig, cost_pf_norm_orig))
    pf_norm_instrdc = np.column_stack((science_pf_norm_instrdc, cost_pf_norm_instrdc))
    pf_norm_instrorb = np.column_stack((science_pf_norm_instrorb, cost_pf_norm_instrorb))
    pf_norm_interinstr = np.column_stack((science_pf_norm_interinstr, cost_pf_norm_interinstr))
    pf_norm_packeff = np.column_stack((science_pf_norm_packeff, cost_pf_norm_packeff))
    pf_norm_spmass = np.column_stack((science_pf_norm_spmass, cost_pf_norm_spmass))
    pf_norm_instrsyn = np.column_stack((science_pf_norm_instrsyn, cost_pf_norm_instrsyn))
    
    hv_orig = compute_hv(pf_norm_orig)
    hv_instrdc = compute_hv(pf_norm_instrdc)
    hv_instrorb = compute_hv(pf_norm_instrorb)
    hv_interinstr = compute_hv(pf_norm_interinstr)
    hv_packeff = compute_hv(pf_norm_packeff)
    hv_spmass = compute_hv(pf_norm_spmass)
    hv_instrsyn = compute_hv(pf_norm_instrsyn)
    
    I_instrdc[i] = hv_instrdc - hv_orig
    I_instrorb[i] = hv_instrorb - hv_orig
    I_interinstr[i] = hv_interinstr - hv_orig
    I_packeff[i] = hv_packeff - hv_orig
    I_spmass[i] = hv_spmass - hv_orig
    I_instrsyn[i] = hv_instrsyn - hv_orig
    
I_op_instrdc = np.mean(I_instrdc)
I_op_instrorb = np.mean(I_instrorb)
I_op_interinstr = np.mean(I_interinstr)
I_op_packeff = np.mean(I_packeff)
I_op_spmass = np.mean(I_spmass)
I_op_instrsyn = np.mean(I_instrsyn)

print('Operator Index - Instrdc : ' + str(I_op_instrdc) + ' +\- ' + str(np.std(I_instrdc)))
print('Operator Index - Instrorb : ' + str(I_op_instrorb) + ' +\- ' + str(np.std(I_instrorb)))
print('Operator Index - Interinstr : ' + str(I_op_interinstr) + ' +\- ' + str(np.std(I_interinstr)))
print('Operator Index - Packeff : ' + str(I_op_packeff) + ' +\- ' + str(np.std(I_packeff)))
print('Operator Index - Spmass : ' + str(I_op_spmass) + ' +\- ' + str(np.std(I_spmass)))
print('Operator Index - Instrsyn : ' + str(I_op_instrsyn) + ' +\- ' + str(np.std(I_instrsyn)))