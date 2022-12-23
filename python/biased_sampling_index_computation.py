# -*- coding: utf-8 -*-
"""
Biased Sampling Impact Index Computation

@author: roshan94
"""
from pygmo import hypervolume
import csv
import numpy as np
#import operator as op
from IPython.core.debugger import set_trace

### Useful functions
def get_csv_filepath_satellite(bias_init_enforced, assigning, run_number):
    # bias_init_enforced = [instrdc, instrorb, interinstr, packeff, spmass, instrsyn, instrcount] boolean array for Assigning problem
    # bias_init_enforced = [instrdc, instrorb, interinstr, packeff, spmass, instrsyn] boolean array for Partitioning problem
    #
    # assigning = True if assigning problem data is to be read, False if partitioning problem data is to be read
    
    filepath = 'C:\\SEAK Lab\\SEAK Lab Github\\VASSAR\\VASSAR_exec_heur\\results\\' # for workstation
    #filepath = 'C:\\Users\\rosha\\Documents\\SEAK Lab Github\\VASSAR\\VASSAR_exec_heur\\results\\' # for laptop
    heurs_list = ['Instrdc','Instrorb','Interinstr','Packeff','Spmass','Instrsyn','Instrcount']
    heur_abbrvs_list = ['d','o','i','p','m','s','c']
    
    filename = 'EpsilonMOEA_emoea_'    
        
    if assigning:
        filepath_prob = 'Assigning\\'
        filename_prob = '_assigning_fullpop.csv'
    else:
        filepath_prob = 'Partitioning\\'
        filename_prob = '_partitioning_fullpop.csv'
        
    filepath2 = ''
    filename2 = ''
    filepath_moea = ''
    filepath_cred = ''
    constr_count = 0
    constraints = 'Bias Init - '
    constraints_abbrv = ''
    for i in range(len(bias_init_enforced)):
        
        if bias_init_enforced[i]:
            constraints = constraints + heurs_list[i]
            constraints_abbrv = constraints_abbrv + heur_abbrvs_list[i]
        else:
            constr_count += 1
        
    if constr_count < len(bias_init_enforced):
        filepath2 = filepath2 + constraints + '\\'
        filename2 = filename2 + constraints_abbrv + 'con2_'
    else:
        filepath_moea = 'Epsilon MOEA\\'
            
    filepath_initialization = ''
    if not assigning:
        filepath_initialization = 'injected initialization\\'
        
    return filepath + filepath_prob + filepath2 + filepath_moea + filepath_cred + filepath_initialization +  filename + str(run_number) + filename2 + filename_prob

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

def extract_data_from_csv(csv_filepath, assigning):
    # intpen_constr_heur = [intpen_constr_instrdc, intpen_constr_instrorb, intpen_constr_interinstr, intpen_constr_packeff, intpen_constr_spmass, intpen_constr_instrsyn[, intpen_constr_instrcount for assigning]] boolean array
    with open(csv_filepath,newline='') as csvfile:
        data = [row for row in csv.reader(csvfile)]
                
        num_func_evals_dat = np.zeros(len(data)-1)
        designs = []
        science_dat = np.zeros(len(data)-1)
        cost_dat = np.zeros(len(data)-1)
        
        valid_count = 0
        for x in range(len(data)-1):
            data_float = list(map(float,data[x+1][1:]))
            if (any(np.isnan(np.array(data_float))) or any(np.isinf(np.array(data_float)))):
                continue
            
            designs.append(data[x+1][0])
            num_func_evals_dat[valid_count] = int(data[x+1][1])
            science_dat[valid_count] = -float(data[x+1][2]) 
            cost_dat[valid_count] = float(data[x+1][3])
            
            valid_count += 1
            
    #archs = archs_dat[:valid_count]
    num_func_evals = num_func_evals_dat[:valid_count]
    science = science_dat[:valid_count]
    cost = cost_dat[:valid_count]
    
    ## Sort num_fun_evals (and objectives and heuristic scores) in ascending order
    n_func_evals = num_func_evals
    sort_indices = np.argsort(n_func_evals)
    #designs_sorted = list(np.array(designs)[sort_indices])
    science_sorted = list(science[sort_indices])
    cost_sorted = list(cost[sort_indices])
    
    nfe_list_sorted = list(n_func_evals[sort_indices])
    
    ## Extract the objectives of the initial population
    
    obj_weights = [0.4, 7250]
    if assigning:
        obj_weights = [0.425, 2.5e4] # normalization weights for objectives
        
    nfe_index_current = find_last_index(0, nfe_list_sorted)
    science_init = science_sorted[:nfe_index_current]
    cost_init = cost_sorted[:nfe_index_current]
    
    init_population = np.zeros((len(science_init), 2))
    
    for i in range(len(science_init)):
        init_population[i,0] = science_init[i]*obj_weights[0]
        init_population[i,1] = cost_init[i]*obj_weights[1]
        
    return init_population

def find_norm_vals(pop):
    science_max = np.amax([x[0] for x in pop])
    science_min = np.amin([x[0] for x in pop])
    cost_max = np.amax([x[1] for x in pop])
    cost_min = np.amin([x[1] for x in pop])
    
    return science_max, science_min, cost_max, cost_min

#### OPERATION
assign_prob = True

if assign_prob:
    eps_moea_enforced = [False, False, False, False, False, False, False]
    instrcount_bias_init_enforced = [False, False, False, False, False, False, True]
else: # Partitioning problem, not used currently
    eps_moea_enforced = [False, False, False, False, False, False]
    instrdc_bias_init_enforced = [True, False, False, False, False, False]
    instrorb_bias_init_enforced = [False, True, False, False, False, False]
    interinstr_bias_init_enforced = [False, False, True, False, False, False]
    packeff_bias_init_enforced = [False, False, False, True, False, False]
    spmass_bias_init_enforced = [False, False, False, False, True, False]
    instrsyn_bias_init_enforced = [False, False, False, False, False, True]

bias_sampling_cases = np.vstack((eps_moea_enforced, instrcount_bias_init_enforced))

n_runs = 10
n_cases = 2

## Extract and store initial populations
init_pop_allcases = {}
for i in range(n_cases):
    init_pop_allruns = {}
    for j in range(n_runs):
        filepath_run = get_csv_filepath_satellite(bias_sampling_cases[i], assign_prob, j)
        init_pop_run = extract_data_from_csv(filepath_run, assign_prob)
        init_pop_allruns['run'+str(j)] = init_pop_run
    init_pop_allcases['case'+str(i)] = init_pop_allruns
        
## Find overall normalization constants
science_max_all = np.zeros((n_cases, n_runs))
science_min_all = np.zeros((n_cases, n_runs))
cost_max_all = np.zeros((n_cases, n_runs))
cost_min_all = np.zeros((n_cases, n_runs))

for i in range(n_cases):
    init_pop_caseruns = init_pop_allcases['case'+str(i)]
    for j in range(n_runs):
        init_pop_run = init_pop_caseruns['run'+str(j)]
        sc_max, sc_min, cost_max, cost_min = find_norm_vals(init_pop_run)
        science_max_all[i][j] = sc_max
        science_min_all[i][j] = sc_min
        cost_max_all[i][j] = cost_max
        cost_min_all[i][j] = cost_min
    
science_norm_max = np.amax(science_max_all)
science_norm_min = np.amin(science_min_all)
cost_norm_max = np.amax(cost_max_all)
cost_norm_min = np.amin(cost_min_all)

#science_norm_max = np.zeros((n_runs))
#science_norm_min = np.zeros((n_runs))
#cost_norm_max = np.zeros((n_runs))
#cost_norm_min = np.zeros((n_runs))

#for i in range(n_runs):
    #science_norm_max[i] = np.amax(science_max_all[:,i])
    #science_norm_min[i] = np.amin(science_min_all[:,i])
    #cost_norm_max[i] = np.amax(cost_max_all[:,i])
    #cost_norm_min[i] = np.amin(cost_min_all[:,i])

## Compute hypervolumes
hv_cases = np.zeros((n_runs, n_cases))
for i in range(n_cases):
    init_pop_caseruns = init_pop_allcases['case'+str(i)]
    for j in range(n_runs):
        init_pop_run = init_pop_caseruns['run'+str(j)]
        science_run = [x[0] for x in init_pop_run]
        cost_run = [x[1] for x in init_pop_run]
        
        science_norm_run = (science_run - science_norm_min)/(science_norm_max - science_norm_min)
        cost_norm_run = (cost_run - cost_norm_min)/(cost_norm_max - cost_norm_min)
        
        #science_norm_run = (science_run - science_norm_min[j])/(science_norm_max[j] - science_norm_min[j])
        #cost_norm_run = (cost_run - cost_norm_min[j])/(cost_norm_max[j] - cost_norm_min[j])
        
        pop_norm = np.column_stack((science_norm_run, cost_norm_run))
        hv_cases[j,i] = compute_hv(pop_norm)
        
## Compute indices
indices = np.zeros((n_runs, n_cases-1))
for i in range(n_cases-1):
    indices = hv_cases[:,i+1] - hv_cases[:,0]
    
print('Instrument Count Impact Index: ' + str(np.mean(indices)) + ' +/- ' + str(np.std(indices)))

### Computing minumum percentile for positive indices 
percentile_vals = np.linspace(1, 100, 100)
if hv_cases.shape[1] == 2:
    n_perctile_heurs = [0]
else:    
    n_perctile_heurs = np.zeros((indices.shape[1]))
    
for i in range(len(n_perctile_heurs)):
    if hv_cases.shape[1] == 2:
        I_current_heur = indices
    else:
        I_current_heur = indices[:,i]
    for j in range(len(percentile_vals)):
        pctile = np.percentile(I_current_heur, percentile_vals[j], method='interpolated_inverted_cdf')
        if pctile > 0:
            n_perctile_heurs[i] = percentile_vals[j]
            break
        if j == len(percentile_vals)-1:
            n_perctile_heurs[i] = percentile_vals[j]
        
