# -*- coding: utf-8 -*-
"""
Hypervolume computation for both satellite problems

@author: roshan94
"""
from pygmo import hypervolume
import csv
import statistics
import numpy as np
#import operator as op
from scipy.stats import mannwhitneyu
from scipy.stats import ttest_ind
from functools import reduce
from itertools import combinations
import matplotlib.pyplot as plt
from IPython.core.debugger import set_trace

assigning_problem = False
num_runs = 30 # number of runs for each case
#threshold_hv = 0.85

credit_assignment = 1 # 0 -> offspring parent dominance, 1 -> set improvement dominance, 2 -> set contribution dominance

#### Useful functions and parameter defintions 
def get_objectives(obj1_array, obj2_array, index):
    return obj1_array[index], obj2_array[index]

def find_last_index(val,search_list):
    if val in search_list:
        idx = len(search_list) - search_list[::-1].index(val) - 1
    else:
        idx = 0
    return idx

def find_closest_index(val,search_list):
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

def get_array_element(array, index):
    return array[index]

#### Determine csv filepath from given case type for one of the satellite problems
def get_csv_filepath_satellite(instrdc_constrained, instrorb_constrained, interinstr_constrained, packeff_constrained, spmass_constrained, instrsyn_constrained, instrcount_constrained, cred_strat, assigning, run_number):
    # instrdc_constrained = [int_pen, AOS, bias_init, ACH] boolean array
    # instrorb_constrained = [int_pen, AOS, bias_init, ACH] boolean array
    # interinstr_constrained = [int_pen, AOS, bias_init, ACH] boolean array
    # packeff_constrained = [int_pen, AOS, bias_init, ACH] boolean array
    # spmass_constrained = [int_pen, AOS, bias_init, ACH] boolean array
    # instrsyn_constrained = [int_pen, AOS, bias_init, ACH] boolean array
    # instrcount_constrained = [int_pen, AOS, bias_init, ACH] boolean array # only for assigning problem
    #
    # assigning = True if assigning problem data is to be read, False if partitioning problem data is to be read
    
    filepath = 'C:\\SEAK Lab\\SEAK Lab Github\\VASSAR\\VASSAR_exec_heur\\results\\' # for workstation
    #filepath = 'C:\\Users\\rosha\\Documents\\SEAK Lab Github\\VASSAR\\VASSAR_exec_heur\\results\\' # for laptop
    methods = ['Int Pen','AOS','Bias Init','ACH']
    heurs_list = ['Instrdc','Instrorb','Interinstr','Packeff','Spmass','Instrsyn','Instrcount']
    heur_abbrvs_list = ['d','o','i','p','m','s','c']
    if assigning:
        heur_bools = np.vstack((instrdc_constrained, instrorb_constrained, interinstr_constrained, packeff_constrained, spmass_constrained, instrsyn_constrained, instrcount_constrained))
    else:
        heur_bools = np.vstack((instrdc_constrained, instrorb_constrained, interinstr_constrained, packeff_constrained, spmass_constrained, instrsyn_constrained))
    aos_bools = [x[1] for x in heur_bools]
    
    if (any(aos_bools)):
        filename = 'AOSMOEA_emoea_'
    else:
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
    for i in range(len(heur_bools[0])):
        constraints = methods[i] + ' - '
        constraints_abbrv = ''
        heur_count = 0
        for j in range(len(heur_bools)):
            if heur_bools[j][i]:
                constraints = constraints + heurs_list[j]
                constraints_abbrv = constraints_abbrv + heur_abbrvs_list[j]
            else:
                heur_count += 1
        
        if heur_count < len(heur_bools):
            filepath2 = filepath2 + constraints + '\\'
            filename2 = filename2 + constraints_abbrv + 'con' + str(i) + '_'
        else:
            constr_count += 1
        
    if (constr_count == len(heur_bools[0])):
        filepath_moea = 'Epsilon MOEA\\'
    else:
        if (any(aos_bools)):
            filepath_cred = 'offspring parent dominance\\'
            if cred_strat == 1:
                filepath_cred = 'set improvement dominance\\'
            elif cred_strat == 2:
                filepath_cred = 'set contribution dominance\\'
            
    filepath_initialization = ''
    if not assigning:
        filepath_initialization = 'injected initialization\\'
        
    return filepath + filepath_prob + filepath2 + filepath_moea + filepath_cred + filepath_initialization +  filename + str(run_number) + filename2 + filename_prob

#### Define NFE array for hypervolume computation (based on number of evaluations in optimization runs)
np.set_printoptions(threshold=np.inf)

# Create array of NFE values at which to compute hypervolume 
# Assumes max function evaluations for a) assigning problem is 5000 and b) partitioning problem is 1000
nfe_increment = 50
    
n_iter_total = 80 # Total number of points in NFE array (1 more than input value to incorporate 0)
n_iter_init = 60 # Number of initial points in NFE array separated by nfe_increment (the rest after that are separated by 2*nfe_increment)
nfe_array = np.zeros(n_iter_total+1)
for i in range(n_iter_init):
    nfe_array[i] = nfe_increment*i
    
for i in range(n_iter_total - n_iter_init + 1):
    nfe_array[n_iter_init+i] = nfe_increment*n_iter_init + (2*nfe_increment)*i
    
#### Extract Pareto Front and normalization constants data from csv file
def extract_data_from_csv(csv_filepath, assigning, intpen_constr_heur):
    # intpen_constr_heur = [intpen_constr_instrdc, intpen_constr_instrorb, intpen_constr_interinstr, intpen_constr_packeff, intpen_constr_spmass, intpen_constr_instrsyn[, intpen_constr_instrcount for assigning]] boolean array
    with open(csv_filepath,newline='') as csvfile:
        data = [row for row in csv.reader(csvfile)]
                
        num_func_evals_dat = np.zeros(len(data)-1)
        designs = []
        science_dat = np.zeros(len(data)-1)
        cost_dat = np.zeros(len(data)-1)
        
        instrdc_scores_dat = np.zeros(len(data)-1)
        instrorb_scores_dat = np.zeros(len(data)-1)
        interinstr_scores_dat = np.zeros(len(data)-1)
        packeff_scores_dat = np.zeros(len(data)-1)
        spmass_scores_dat = np.zeros(len(data)-1)
        instrsyn_scores_dat = np.zeros(len(data)-1)
        if assigning:
            instrcount_scores_dat = np.zeros(len(data)-1)
        
        valid_count = 0
        for x in range(len(data)-1):
            data_float = list(map(float,data[x+1][1:]))
            if (any(np.isnan(np.array(data_float))) or any(np.isinf(np.array(data_float)))):
                continue
            
            designs.append(data[x+1][0])
            num_func_evals_dat[valid_count] = int(data[x+1][1])
            science_dat[valid_count] = -float(data[x+1][2]) 
            cost_dat[valid_count] = float(data[x+1][3])
            
            instrdc_scores_dat[valid_count] = float(data[x+1][4])
            instrorb_scores_dat[valid_count] = float(data[x+1][5])
            interinstr_scores_dat[valid_count] = float(data[x+1][6])
            packeff_scores_dat[valid_count] = float(data[x+1][7])
            spmass_scores_dat[valid_count] = float(data[x+1][8])
            instrsyn_scores_dat[valid_count] = float(data[x+1][9])
            if assigning:
                instrcount_scores_dat[valid_count] = float(data[x+1][10])
    
            valid_count += 1
            
    #archs = archs_dat[:valid_count]
    num_func_evals = num_func_evals_dat[:valid_count]
    science = science_dat[:valid_count]
    cost = cost_dat[:valid_count]
    
    instrdc_scores = instrdc_scores_dat[:valid_count]
    instrorb_scores = instrorb_scores_dat[:valid_count]
    interinstr_scores = interinstr_scores_dat[:valid_count]
    packeff_scores = packeff_scores_dat[:valid_count]
    spmass_scores = spmass_scores_dat[:valid_count]
    instrsyn_scores = instrsyn_scores_dat[:valid_count]
    if assigning:
        instrcount_scores = instrcount_scores_dat[:valid_count]
            
    #print('science')
    #print(science)
    #print('\n')
    #print('cost')
    #print(cost)
    #print('\n')
    
    ## Sort num_fun_evals (and objectives and heuristic scores) in ascending order
    n_func_evals = num_func_evals
    sort_indices = np.argsort(n_func_evals)
    designs_sorted = list(np.array(designs)[sort_indices])
    science_sorted = list(science[sort_indices])
    cost_sorted = list(cost[sort_indices])
    
    instrdc_scores_sorted = list(instrdc_scores[sort_indices])
    instrorb_scores_sorted = list(instrorb_scores[sort_indices])
    interinstr_scores_sorted = list(interinstr_scores[sort_indices])
    packeff_scores_sorted = list(packeff_scores[sort_indices])
    spmass_scores_sorted = list(spmass_scores[sort_indices])
    instrsyn_scores_sorted = list(instrsyn_scores[sort_indices])
    if assigning:
        instrcount_scores_sorted = list(instrcount_scores[sort_indices])
    
    #archs_sorted = []
    #for i in range(len(sort_indices)):
        #archs_sorted.append(archs[sort_indices[i]])
    
    if assigning:
        heur_scores_sorted = np.vstack((instrdc_scores_sorted, instrorb_scores_sorted, interinstr_scores_sorted, packeff_scores_sorted, spmass_scores_sorted, instrsyn_scores_sorted, instrcount_scores_sorted))
    else:
        heur_scores_sorted = np.vstack((instrdc_scores_sorted, instrorb_scores_sorted, interinstr_scores_sorted, packeff_scores_sorted, spmass_scores_sorted, instrsyn_scores_sorted))
    
    ## Compute objectives (specifically for the interior penalty cases)
    if assigning:
        heur_objs_norm = [0.425, 2.5e4] # normalization weights for objectives
    else:
        heur_objs_norm = [0.4, 7250]
        
    heur_weight = 1 # change to 0.1 for interior penalty for either problem
    
    heur_objs = np.zeros(len(instrdc_scores_sorted))
    if any(intpen_constr_heur):
        heur_index_array = np.arange(len(intpen_constr_heur))
        heur_index_used = [v for i, v in enumerate(heur_index_array) if intpen_constr_heur[i] == True]
     
        for idx in heur_index_used:
            heur_scores_idx = heur_scores_sorted[idx]
            heur_objs_idx = [heur_weight*x for x in heur_scores_idx]
            heur_objs = np.add(heur_objs, heur_objs_idx)
            
        heur_objs = np.divide(heur_objs, len(heur_index_used))
    
    science_objs_sorted = list(np.subtract(science_sorted, heur_objs))
    cost_objs_sorted = list(np.subtract(cost_sorted, heur_objs))
    
    true_science = [k*heur_objs_norm[0] for k in science_objs_sorted] 
    true_cost = [k*heur_objs_norm[1] for k in cost_objs_sorted]
    
    ## Determine normalizing objective scores and compute pareto fronts for penalized and true objectives as well as for true objectives of only feasible designs 
    nfe_list_sorted = list(n_func_evals[sort_indices])
    
    #max_func_evals = nfe_list_sorted[-1]
    max_func_evals = 5000 # some runs for some cases run upto 5001 evaluations, which causes hv array length issues

    pareto_front_dict = {}
    #pareto_front_instrdc_dict = {}
    #pareto_front_instrorb_dict = {}
    #pareto_front_interinstr_dict = {}
    #pareto_front_packeff_dict = {}
    #pareto_front_spmass_dict = {}
    #pareto_front_instrsyn_dict = {}
    #pareto_front_archs_dict = {}

    pf_normalize_max_obj1 = []
    pf_normalize_min_obj1 = []
    pf_normalize_max_obj2 = []
    pf_normalize_min_obj2 = []
    
    count = 0
    pop_size = int(find_last_index(0, nfe_list_sorted))

    for i in range(len(nfe_array)):
        #print('iter = ' + str(i))
        nfe_val = nfe_array[i]
    
        if (nfe_list_sorted[0] == 0):
            if (nfe_val <= pop_size): # population size set at 300 in java code, but maybe different here due to NaNs
                nfe_index_current = pop_size
                nfe_index_previous = 0
            else:
                nfe_index_current = find_closest_index(nfe_val, nfe_list_sorted)
                nfe_index_previous = find_closest_index(nfe_array[i-1], nfe_list_sorted)   
        else:
            if (nfe_val <= nfe_list_sorted[0]):
                nfe_index_previous = 0 
                nfe_index_current = find_closest_index(0, nfe_list_sorted)
            else:
                nfe_index_current = find_closest_index(nfe_val, nfe_list_sorted)
                nfe_index_previous = find_closest_index(nfe_array[i-1], nfe_list_sorted)
                
        if (nfe_index_current == nfe_index_previous):
            nfe_index_current = nfe_index_current + 1
        
        nfe_array_current = nfe_list_sorted[nfe_index_previous:nfe_index_current]
        current_population = []
        for j in range(len(nfe_array_current)):
            current_population.append([true_science[nfe_index_previous+j], true_cost[nfe_index_previous+j]])
            
        #if ("AOS - Orient\\" in csv_filepath) and ("emoea_16" in csv_filepath):
            #set_trace()

        if (i != 0):
            previous_pareto_front = pareto_front_dict[nfe_array[i-1]].tolist()
            for k2 in range(len(previous_pareto_front)):
                current_population.append(previous_pareto_front[k2])
        
        current_pareto_front_all = compute_pareto_front(current_population)
        #current_pareto_front = list(set(current_pareto_front_all))
        current_pareto_front = np.unique(current_pareto_front_all, axis=0)
    
        #current_pareto_instrdc_scores = []
        #current_pareto_instrorb_scores = []
        #current_pareto_interinstr_scores = []
        #current_pareto_packeff_scores = []
        #current_pareto_spmass_scores = []
        #current_pareto_instrsyn_scores = []
        #current_pareto_archs = []
        ##current_pareto = []
        ##for pareto_arch in current_pareto_front:
            ##arch_index = true_science.index(pareto_arch[0])
            #arch_instrdc_score = get_array_element(instrdc_scores_sorted, arch_index)
            #arch_instrorb_score = get_array_element(instrorb_scores_sorted, arch_index)
            #arch_interinstr_score = get_array_element(interinstr_scores_sorted, arch_index)
            #arch_packeff_score = get_array_element(packeff_scores_sorted, arch_index)
            #arch_spmass_score = get_array_element(spmass_scores_sorted, arch_index)
            #arch_instrsyn_score = get_array_element(instsyn_scores_sorted, arch_index)
            
            #current_pareto_instrdc_scores.append(arch_instrdc_score)
            #current_pareto_instrorb_scores.append(arch_instrorb_score)
            #current_pareto_interinstr_scores.append(arch_interinstr_score)
            #current_pareto_packeff_scores.append(arch_packeff_score)
            #current_pareto_spmass_scores.append(arch_spmass_score)
            #current_pareto_instrsyn_scores.append(arch_instrsyn_score)
            #current_pareto_archs.append(get_array_element(archs_sorted, arch_index))
            
            ##science_arch, cost_arch = get_objectives(science_sorted, cost_sorted, arch_index)
            #set_trace()
            ##current_pareto.append([science_arch, cost_arch])
                   
        pareto_front_dict[nfe_val] = current_pareto_front
        #pareto_front_instrdc_dict[nfe_val] = current_pareto_instrdc_scores
        #pareto_front_instrorb_dict[nfe_val] = current_pareto_instrorb_scores
        #pareto_front_interinstr_dict[nfe_val] = current_pareto_interinstr_scores
        #pareto_front_packeff_dict[nfe_val] = current_pareto_packeff_scores
        #pareto_front_spmass_dict[nfe_val] = current_pareto_spmass_scores
        #pareto_front_instrsyn_dict[nfe_val] = current_pareto_instrsyn_scores
        #pareto_front_archs_dict[nfe_val] = 

        pf_nfeval_obj1 = [row[0] for row in current_pareto_front]
        pf_nfeval_obj2 = [row[1] for row in current_pareto_front]
        
        #pf_objs_nfeval_obj1 = [row[0] for row in current_pareto_objs]
        #pf_objs_nfeval_obj2 = [row[1] for row in current_pareto_objs]
        
        pf_normalize_max_obj1.append(np.max(pf_nfeval_obj1))
        pf_normalize_min_obj1.append(np.min(pf_nfeval_obj1))
        pf_normalize_max_obj2.append(np.max(pf_nfeval_obj2))
        pf_normalize_min_obj2.append(np.min(pf_nfeval_obj2))
        
        #pf_objs_normalize_max_obj1.append(np.max(pf_objs_nfeval_obj1))
        #pf_objs_normalize_min_obj1.append(np.min(pf_objs_nfeval_obj1))
        #pf_objs_normalize_max_obj2.append(np.max(pf_objs_nfeval_obj2))
        #pf_objs_normalize_min_obj2.append(np.min(pf_objs_nfeval_obj2))

    ### Computing obj_normalize_fullrun, obj_normalize_afterjump and obj_normalized_true_fullrun using the entire run 
    #obj_normalize_max_fullrun = [np.max(pen_obj1_constr_sorted), np.max(pen_obj2_constr_sorted)]
    #obj_normalize_min_fullrun = [np.min(pen_obj1_constr_sorted), np.min(pen_obj2_constr_sorted)]
    
    #obj_true_normalize_max_fullrun = [np.max(true_obj1_sorted), np.max(true_obj2_sorted)]
    #obj_true_normalize_min_fullrun = [np.min(true_obj1_sorted), np.min(true_obj2_sorted)]

    #obj_normalize_fullrun = [obj_normalize_min_fullrun, obj_normalize_max_fullrun]
    
    #obj_normalize_true_fullrun = [obj_true_normalize_min_fullrun, obj_true_normalize_max_fullrun]
    
    ### Computing obj_normalize_fullrun using the pareto fronts
    obj_normalize_max_fullrun = [np.max(pf_normalize_max_obj1), np.max(pf_normalize_max_obj2)]
    obj_normalize_min_fullrun = [np.min(pf_normalize_min_obj1), np.min(pf_normalize_min_obj2)]  
    
    #obj_objs_normalize_max_fullrun = [np.max(pf_objs_normalize_max_obj1), np.max(pf_objs_normalize_max_obj2)]
    #obj_objs_normalize_min_fullrun = [np.min(pf_objs_normalize_min_obj1), np.min(pf_objs_normalize_min_obj2)]
    
    obj_normalize_fullrun = [obj_normalize_min_fullrun, obj_normalize_max_fullrun]
    
    #obj_normalize_objs_fullrun = [obj_objs_normalize_min_fullrun, obj_objs_normalize_max_fullrun]
    
    #return pareto_front_dict, pareto_front_obj_dict, obj_normalize_fullrun, obj_normalize_objs_fullrun, max_func_evals
    return pareto_front_dict, obj_normalize_fullrun, max_func_evals

#### Compute overall normalization objectives for single case study/all compared case studies (from the complete runs)
def compute_overall_norm_objs(objs_normalization_objs):
    # Each input is a dictionary with key as the case study/run string and value as the corresponding 2D array
    
    obj1_max_objs_allcases = np.zeros(len(objs_normalization_objs))
    obj1_min_objs_allcases = np.zeros(len(objs_normalization_objs))
    obj2_max_objs_allcases = np.zeros(len(objs_normalization_objs))
    obj2_min_objs_allcases = np.zeros(len(objs_normalization_objs))
        
    i = 0
    for key in objs_normalization_objs:
        current_objs_norm_objs = objs_normalization_objs[key]
        
        obj1_max_objs_allcases[i] = current_objs_norm_objs[1][0]
        obj2_max_objs_allcases[i] = current_objs_norm_objs[1][1]
        obj1_min_objs_allcases[i] = current_objs_norm_objs[0][0]
        obj2_min_objs_allcases[i] = current_objs_norm_objs[0][1]
        i += 1
    
    obj1_min_objs_overall = np.min(obj1_min_objs_allcases)
    obj2_min_objs_overall = np.min(obj2_min_objs_allcases)
    obj1_max_objs_overall = np.max(obj1_max_objs_allcases)
    obj2_max_objs_overall = np.max(obj2_max_objs_allcases)
            
    obj_norm_objs_overall = [[obj1_min_objs_overall, obj2_min_objs_overall], [obj1_max_objs_overall, obj2_max_objs_overall]]    
     
    return obj_norm_objs_overall

#### Compute hypervolume arrays from copmuted pareto fronts and normalization constants
def compute_hv_arrays_from_csv_data(pf_dict, obj_norm_full, max_fun_evals):
    obj_norm_min_full = obj_norm_full[0]
    obj_norm_max_full = obj_norm_full[1]

    ## Normalize the pareto front objectives and compute the hypervolume
    hypervol_full_dict = []

    for nfe_val in nfe_array:
        #print('iter = ' + str(nfe_val))
    
        current_pareto_front = pf_dict[nfe_val]
        current_pf_normalized = []
        current_pf_objs_normalized = []
        for pareto_design in current_pareto_front:
            obj1_normalized = (pareto_design[0] - obj_norm_min_full[0])/(obj_norm_max_full[0] - obj_norm_min_full[0])
            obj2_normalized = (pareto_design[1] - obj_norm_min_full[1])/(obj_norm_max_full[1] - obj_norm_min_full[1])
            current_pf_normalized.append([obj1_normalized, obj2_normalized])
            
        current_hv = compute_hv(current_pf_normalized)
        hypervol_full_dict.append([nfe_val, current_hv])
        
    return hypervol_full_dict

### Compute the max hypervolume of all cases
def get_max_hv(hypervol_dict_allcases):
    hv_max_allcases = np.zeros((len(hypervol_dict_allcases), len(hypervol_dict_allcases['case'+str(0)])))
    for i in range(len(hypervol_dict_allcases)):
        hv_dict_allruns_case = hypervol_dict_allcases['case'+str(i)]
        for j in range(len(hv_dict_allruns_case)):
            hv_vals_run = [x[1] for x in hv_dict_allruns_case['run'+str(j)]]
            hv_max_run = np.amax(hv_vals_run)
            hv_max_allcases[i,j] = hv_max_run
        
    return np.amax(hv_max_allcases)

### Compute array of NFE values for reaching threshold hypervolume for different runs of a particular case
def compute_nfe_hypervolume_attained(hv_dict, hv_threshold):
    #hv_threshold = 0.75 # Threshold HV value to reach, user parameter
    n_runs = len(hv_dict)
    nfe_hv_attained = []
    for key in hv_dict:
        hv_array_run = hv_dict[key]
        nfe_array_run = [hv_array[0] for hv_array in hv_array_run]
        hv_val_array = [hv_array[1] for hv_array in hv_array_run]
        nfe_hv_attained_run = nfe_array_run[-1] + 100
        for i in range(len(hv_val_array)):
            if (hv_val_array[i] >= hv_threshold):
                nfe_hv_attained_run = nfe_array_run[i]
                break
                
        #hv_val_diff = [np.abs(x - hv_threshold) for x in hv_val_array]
        #index = np.argmin(hv_val_diff)
        nfe_hv_attained.append(nfe_hv_attained_run)
        
    return nfe_hv_attained
    
### Plot fraction of runs attaining threshold hypervolume vs NFE
def plot_fraction_hypervolume_attained(nfe_hv_attained_dict, nfe_array, colour_array, linestyles_array, casename_array, savefig_name):
    fig1 = plt.figure()
    n_cases = len(nfe_hv_attained_dict)
    case_idx = 0
    for case_key in nfe_hv_attained_dict:
        nfe_hv_attained_case = nfe_hv_attained_dict[case_key]
        n_runs = len(nfe_hv_attained_case)
        frac_runs_hv_attained = np.zeros(len(nfe_array))
        for i in range(len(nfe_array)):
            idx_runs_hv_attained = [idx for idx, val in enumerate(nfe_hv_attained_case) if val <= nfe_array[i]]
            frac_runs_hv_attained[i] = len(idx_runs_hv_attained)/n_runs
            
        plt.plot(nfe_array, frac_runs_hv_attained, linestyle=linestyles_array[case_idx], color=colour_array[case_idx], label=casename_array[case_idx])
        case_idx += 1
    
    plt.xlabel(r'Number of Function Evaluations',fontsize=14)
    plt.ylabel(r'Fraction of runs HV $\geq$ 0.8 $\times$ max HV',fontsize=14)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.legend(loc='upper center', bbox_to_anchor=(0.5,1.10), ncol=3, borderaxespad=0 ,prop={"size":12})
    plt.show()
    fig1.savefig('frac_hv_attained_' + savefig_name + '.pdf', format='pdf')
    
def compute_hypervolume_stats(hypervols_dict):
    hv_dict_keys = list(hypervols_dict.keys())
    hv_dict_0 = hypervols_dict[hv_dict_keys[0]]
    nfe_array_0 = [hv_array[0] for hv_array in hv_dict_0]
    n_datapoints = len(nfe_array_0)
    hypervol_median = np.zeros(n_datapoints)
    hypervol_1q = np.zeros(n_datapoints)
    hypervol_3q = np.zeros(n_datapoints)
    for i in range(n_datapoints):
        hypervol_vals = []
        for key in hypervols_dict:
            hv_dict_j = hypervols_dict[key]
            hv_current_array = [hv_array[1] for hv_array in hv_dict_j]
            hypervol_vals.append(hv_current_array[i])
        hypervol_median[i] = statistics.median(hypervol_vals)
        #hypervol_median[i] = statistics.mean(hypervol_vals)
        hypervol_1q[i] = np.percentile(hypervol_vals, 25)
        #hypervol_1q[i] = hypervol_median[i] - statistics.stdev(hypervol_vals)
        hypervol_3q[i] = np.percentile(hypervol_vals, 75)
        #hypervol_3q[i] = hypervol_median[i] + statistics.stdev(hypervol_vals)
        
    return hypervol_median, hypervol_1q, hypervol_3q, nfe_array_0

def plot_hypervolume_stats(hv_median_case, hv_1q_case, hv_3q_case, nfe_array, savefig_name):
    fig1 = plt.figure(1)
    plt.plot(nfe_array,hv_median_case, 'b-', label='Median')
    plt.plot(nfe_array,hv_1q_case, 'r-', label='1st Quartile')
    plt.plot(nfe_array,hv_3q_case, 'g-', label='3rd Quartile')
    plt.xlabel('Number of Function Evaluations')
    plt.ylabel('Hypervolume')
    plt.title('Averaged Hypervolume vs NFE')
    plt.legend(loc='lower right')
    plt.show()
    #fig1.savefig('HV_plot_averaged_' + savefig_name + '.png')
    
def plot_hypervolume_stats_allcases(hv_median_dict, hv_1q_dict, hv_3q_dict, nfe_array, colour_array, alpha_array, linestyles_array, hatches_case, casename_array, plot_title, savefig_name, incl_legend, use_ylim):
    fig1 = plt.figure()
    number_cases = len(hv_median_dict)
    #print('n_cases')
    #print(number_cases)
    for i in range(number_cases):
        #print(print(marker_array[i]+'*'))
        if all(h is None for h in hatches_case):
            plt.plot(nfe_array, hv_median_dict['case'+str(i)], linewidth=2.5, color=colour_array[i], linestyle=linestyles_array[i], label=casename_array[i])
            plt.fill_between(nfe_array, hv_1q_dict['case'+str(i)], hv_3q_dict['case'+str(i)], color=colour_array[i], hatch=hatches_case[i], linestyle=linestyles_array[i], edgecolor="none", alpha=alpha_array[i])
        else:
            plt.plot(nfe_array, hv_median_dict['case'+str(i)], linewidth=2.5, linestyle=linestyles_array[i], label=casename_array[i])
            plt.fill_between(nfe_array, hv_1q_dict['case'+str(i)], hv_3q_dict['case'+str(i)], facecolor="none", edgecolor=colour_array[i], hatch=hatches_case[i], linestyle=linestyles_array[i], alpha=alpha_array[i])
        
        #plt.plot(nfe_array, hv_1q_dict['case'+str(i)], '--', color=colour_array[i])#, label=casename_array[i]+' 1st Quartile')
        #plt.plot(nfe_array, hv_3q_dict['case'+str(i)], '--', color=colour_array[i])#, label=casename_array[i]+' 3rd Quartile')
    if use_ylim:
        plt.ylim(0.0, 1.0) # reviewer requested modification
    plt.xlabel(r'Number of Function Evaluations',fontsize=13)
    plt.ylabel(r'Hypervolume',fontsize=13)
    plt.xticks(fontsize=13)
    plt.yticks(fontsize=13)
    #plt.title(plot_title)
    if incl_legend:
        plt.legend(loc='upper center', bbox_to_anchor=(0.5,1.10), ncol=3, borderaxespad=0, prop={"size":12})
    plt.show()
    fig1.savefig('HV_plot_averaged_' + savefig_name + '.pdf', format='pdf')
    
def compute_mann_whitney_Uvals(assigning_prob, hv_dict_allcases, nfe_array): # Wilcoxon Rank Sum Test
    # hv_med_array_allcases is a dictionary of length = number of cases
    
    nfe_samples_array = [0, 250, 500, 750, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000]
    n_samples = len(nfe_samples_array)
    
    #n_samples = 11
    #linspace_samples_array = np.linspace(0,1,n_samples)
    #nfe_samples_array = np.floor(np.multiply(linspace_samples_array, nfe_array[-1]))
       
    #n_samples = 20
    #linspace_samples_array = np.linspace(0,1,n_samples)
    #nfe_samples_array = np.floor(np.multiply(linspace_samples_array, nfe_array[-1]))
    nfe_samples_indices_array = np.zeros(len(nfe_samples_array))
    for i in range(len(nfe_samples_array)):
        nfe_samples_indices_array[i] = find_closest_index(nfe_samples_array[i], nfe_array)
        
    hv_samples_allcases_allruns = {}
    for j in range(len(hv_dict_allcases)):
        hv_dict_allruns_currentcase = hv_dict_allcases['case'+str(j)]
        hv_samples_allruns = {}
        for k in range(len(nfe_samples_indices_array)):
            hv_samples_nfe_allruns = np.zeros(len(hv_dict_allruns_currentcase))
            for m in range(len(hv_dict_allruns_currentcase)):
                hv_run = hv_dict_allruns_currentcase['run'+str(m)]
                hv_samples_nfe_allruns[m] = hv_run[int(nfe_samples_indices_array[k])][1]
            hv_samples_allruns['nfe:'+str(int(nfe_samples_array[k]))] = hv_samples_nfe_allruns
        hv_samples_allcases_allruns['case'+str(j)] = hv_samples_allruns
        
    cases_array = np.arange(len(hv_dict_allcases))
    case_combinations = list(combinations(cases_array,2))
    
    U_test_cases = {}
    
    for n in range(len(case_combinations)):
        case_string_x = 'case' + str(case_combinations[n][0])
        case_string_y = 'case' + str(case_combinations[n][1])
        
        hv_allruns_casex = hv_samples_allcases_allruns[case_string_x]
        hv_allruns_casey = hv_samples_allcases_allruns[case_string_y]
        
        U_test_cases_allnfes = {}
        for p in range(len(nfe_samples_indices_array)):
            hv_samples_nfe_casex = hv_allruns_casex['nfe:'+str(int(nfe_samples_array[p]))]
            hv_samples_nfe_casey = hv_allruns_casey['nfe:'+str(int(nfe_samples_array[p]))]
            
            U1, p_val = mannwhitneyu(hv_samples_nfe_casex, hv_samples_nfe_casey, alternative='less')
            t_val, p_val_t = ttest_ind(hv_samples_nfe_casex, hv_samples_nfe_casey, equal_var=False, alternative='less')
            
            U2 = len(hv_samples_nfe_casex)*len(hv_samples_nfe_casey) - U1
            
            U_test = np.min(np.array([U1, U2]))
            
            U_test_cases_allnfes['nfe:'+str(int(nfe_samples_array[p]))] = [U_test, p_val, p_val_t]
        
        dict_key = case_string_x + ' and ' + case_string_y
        U_test_cases[dict_key] = U_test_cases_allnfes
            
    return U_test_cases

#### Define functions to compute and plot hypervolume for single case and all cases
def hypervolume_computation_single_case(case_booleans, prob_assigning, cred_assign, run_nums, case_name):
    ## Computing the pareto fronts and normalization objectives for each run
    obj_norm_allruns = {}
    pf_allruns = {}
    max_f_evals_allruns = np.zeros(run_nums)
    for i in range(run_nums):
        print('Computing Pareto Fronts for run ' + str(i))
        current_csvpath = get_csv_filepath_satellite(case_booleans[:4], case_booleans[4:8], case_booleans[8:12], case_booleans[12:16], case_booleans[16:20], case_booleans[20:24], case_booleans[-1], cred_assign, prob_assigning, i)
        heur_intpen_constr = [case_booleans[0], case_booleans[4], case_booleans[8], case_booleans[12]]
        pf_dict_i, obj_norm_full_i, max_fun_evals_i = extract_data_from_csv(current_csvpath, prob_assigning, heur_intpen_constr, case_booleans[-1])
        pf_allruns['run'+str(i)] = pf_dict_i
        obj_norm_allruns['run'+str(i)] = obj_norm_full_i
        max_f_evals_allruns[i] = max_fun_evals_i
    
    #print('pf_allruns')
    #print(pf_allruns)
    #print('\n')
    ## Use computed normalization objectives and find the overall normalization objectives across all runs
    print('Computing overall normalization constants')
    norm_objs_full_overall = compute_overall_norm_objs(obj_norm_allruns)
    #print('norm_objs_full_overall')
    #print(norm_objs_full_overall)
    #print('\n')
    
    ## Compute Hypervolume values for each run
    hv_dict_allruns = {}
    for j in range(run_nums):
        print('Computing hypervolume values for run ' + str(j))
        hv_dict_j = compute_hv_arrays_from_csv_data(pf_allruns['run'+str(j)], norm_objs_full_overall, max_f_evals_allruns[j])
        hv_dict_allruns['run'+str(j)] = hv_dict_j
        
    #print('hv_dict_allruns')
    #print(hv_dict_allruns)
        
    ## Plotting
    print('Plotting')
    hv_median_all, hv_1q_all, hv_3q_all, nfe_array = compute_hypervolume_stats(hv_dict_allruns)
    plot_hypervolume_stats(hv_median_all, hv_1q_all, hv_3q_all, nfe_array, case_name+'_full')
    
    
def hypervolume_computation_all_cases(case_bools_dict, cred_assign, prob_assigning, run_nums, marker_colours, alpha_vals, case_names):
    num_cases = len(case_bools_dict) # number of cases to compare 

    ## Computing the pareto fronts and normalization objectives for each run in each case
    pf_allcases = {}
    obj_norm_allcasesandruns = {}
    max_f_evals_allcases = {}
    for i in range(num_cases):
        print('Computing Pareto Fronts for runs in Case '+str(i))
        current_case_bools = case_bools_dict['case'+str(i+1)]
        #set_trace()
        pf_allruns_i = {}
        max_f_evals_allruns = np.zeros(run_nums)
        for j in range(run_nums):
            print('Run '+str(j))
            #if np.isin(j, np.array([8])):
                #print('Debug')
            if prob_assigning:
                current_csvpath = get_csv_filepath_satellite(current_case_bools[:4], current_case_bools[4:8], current_case_bools[8:12], current_case_bools[12:16], current_case_bools[16:20], current_case_bools[20:24], current_case_bools[24:28], cred_assign, prob_assigning, j)
            else:
                current_csvpath = get_csv_filepath_satellite(current_case_bools[:4], current_case_bools[4:8], current_case_bools[8:12], current_case_bools[12:16], current_case_bools[16:20], current_case_bools[20:24], current_case_bools[20:24], cred_assign, prob_assigning, j)
            #set_trace()
            if prob_assigning:
                heur_intpen_constr = [current_case_bools[0], current_case_bools[4], current_case_bools[8], current_case_bools[12], current_case_bools[16], current_case_bools[20], current_case_bools[24]]
            else:
                heur_intpen_constr = [current_case_bools[0], current_case_bools[4], current_case_bools[8], current_case_bools[12], current_case_bools[16], current_case_bools[20]]
            pf_dict_j, obj_norm_full_j, max_fun_evals_j = extract_data_from_csv(current_csvpath, prob_assigning, heur_intpen_constr)
            pf_allruns_i['run'+str(j)] = pf_dict_j
            obj_norm_allcasesandruns['case'+str(i+1)+'run'+str(j)] = obj_norm_full_j
            max_f_evals_allruns[j] = max_fun_evals_j
        pf_allcases['case'+str(i+1)] = pf_allruns_i
        max_f_evals_allcases['case'+str(i+1)] = max_f_evals_allruns
    
    ## Use computed normalization objectives and find the overall normalization objectives across all runs and cases
    print('Computing overall normalization constants')
    norm_objs_full_overall = compute_overall_norm_objs(obj_norm_allcasesandruns)
    
    #set_trace()
    
    ## Compute Hypervolume values for each run in each case
    hv_dict_allcases = {}
    hv_dict_median_allcases = {}
    hv_dict_1q_allcases = {}
    hv_dict_3q_allcases = {}
    nfe_array_hv_attained_dict = {}
    for i in range(num_cases):
        print('Computing hypervolume values for runs in Case '+str(i))
        pfs_case_i = pf_allcases['case'+str(i+1)]
        max_func_evals_i = max_f_evals_allcases['case'+str(i+1)]
        hv_dict_allruns = {}
        for j in range(run_nums):
            print('Run '+str(j))
            hv_dict_j = compute_hv_arrays_from_csv_data(pfs_case_i['run'+str(j)], norm_objs_full_overall, max_func_evals_i[j])
            hv_dict_allruns['run'+str(j)] = hv_dict_j
            
        hv_dict_allcases['case'+str(i)] = hv_dict_allruns
        
    # Compute max hypervolume over all cases for nfe_hypervolume_attained
    hv_max = get_max_hv(hv_dict_allcases)
    hv_thresh = 0.8*hv_max
            
    for i in range(num_cases):
        hv_dict_allruns = hv_dict_allcases['case'+str(i)]
        
        print('Computing array of NFE for attaining threshold hypervolume')
        nfe_hv_attained_case = compute_nfe_hypervolume_attained(hv_dict_allruns, hv_thresh)
        nfe_array_hv_attained_dict['case'+str(i)] = nfe_hv_attained_case
                
        print('Computing hypervolume stats')
        hv_med_i, hv_1q_i, hv_3q_i, nfe_array_i = compute_hypervolume_stats(hv_dict_allruns)
                
        hv_dict_median_allcases['case'+str(i)] = hv_med_i
        hv_dict_1q_allcases['case'+str(i)] = hv_1q_i
        hv_dict_3q_allcases['case'+str(i)] = hv_3q_i
        
    print('Computing Wilcoxon Test Statistics')
    U_test_dict = compute_mann_whitney_Uvals(prob_assigning, hv_dict_allcases, nfe_array_i)
          
    return nfe_array_hv_attained_dict, hv_dict_median_allcases, hv_dict_1q_allcases, hv_dict_3q_allcases, U_test_dict, nfe_array_i
    
def plotting_all_cases(nfe_hv_attained_dict, hv_dict_med_allcases, hv_dict_1stq_allcases, hv_dict_3rdq_allcases, nfe_array0, mark_colors, alphas, line_styles_case, hatches_case, names_cases):
    print('Plotting')
    plot_fraction_hypervolume_attained(nfe_hv_attained_dict, nfe_array0, mark_colors, line_styles_case, names_cases, 'allcases')
    plot_hypervolume_stats_allcases(hv_dict_med_allcases, hv_dict_1stq_allcases, hv_dict_3rdq_allcases, nfe_array, mark_colors, alphas, line_styles_case, hatches_case, names_cases, 'Hypervolume', 'allcases', True, True)
    
def plot_upto_nfe(nfe_start, nfe_end, hv_med_allcases, hv_1stq_allcases, hv_3rdq_allcases, nfe_array0, mark_colors, alphas, line_styles_case, hatches_case, names_cases):
    nfe_start_index = find_closest_index(nfe_start, nfe_array0) # find index of closest value in nfe_array to nfe_stert
    nfe_end_index = find_closest_index(nfe_end, nfe_array0) # find index of closest value in nfe_array to nfe_end
    
    nfe_array_lim = nfe_array0[nfe_start_index:nfe_end_index]
    hv_med_lim_allcases = {}
    hv_1q_lim_allcases = {}
    hv_3q_lim_allcases = {}
    for i in range(len(hv_med_allcases)):
        hv_med_lim_allcases['case'+str(i)] = hv_med_allcases['case'+str(i)][nfe_start_index:nfe_end_index]
        hv_1q_lim_allcases['case'+str(i)] = hv_1stq_allcases['case'+str(i)][nfe_start_index:nfe_end_index]
        hv_3q_lim_allcases['case'+str(i)] = hv_3rdq_allcases['case'+str(i)][nfe_start_index:nfe_end_index]
    
    plot_hypervolume_stats_allcases(hv_med_lim_allcases, hv_1q_lim_allcases, hv_3q_lim_allcases, nfe_array_lim, mark_colors, alphas, line_styles_case, hatches_case, names_cases, '', 'hv_upto_nfe', False, False)
    
#######################################################################################################################################################################################################################################################################################################################################################################################################
    # PROGRAM OPERATION
#######################################################################################################################################################################################################################################################################################################################################################################################################

## NOTE: For Int Pen cases, change heur_weight on line 261 to appropriate value

#### Comparing Simple E-MOEA with AOS - all heuristics and AOS - promising heuristics
cases_dict = {}

# bools = [int_pen_instrdc, AOS_instrdc, bias_init_instrdc, ACH_instrdc, int_pen_instrorb, AOS_instrorb, bias_init_instrorb, ACH_instrorb, int_pen_interinstr, AOS_interinstr, bias_init_interinstr, ACH_interinstr, int_pen_packeff, AOS_packeff, bias_init_packeff, ACH_packeff, int_pen_spmass, AOS_spmass, bias_init_spmass, ACH_spmass, int_pen_instrsyn, AOS_instrsyn, bias_init_instrsyn, ACH_instrsyn, int_pen_instrcount, AOS_instrcount, bias_init_instrcount, ACH_instrcount]
case1_bools = [False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False] # Simple E-MOEA
case2_bools = [False, True, False, False, False, True, False, False, False, True, False, False, False, True, False, False, False, True, False, False, False, True, False, False, False, True, False, False] #  AOS - all heuristics
#case2_bools = [True, False, False, False, True, False, False, False, True, False, False, False, True, False, False, False, True, False, False, False, True, False, False, False, True, False, False, False] #  Int Pen - all heuristics
#case2_bools = [False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, True, False] #  Bias Init - Instrcount
if assigning_problem:
    #case3_bools = [False, True, False, False, False, False, False, False, False, True, False, False, False, False, False, False, False, True, False, False, False, False, False, False, False, True, False, False] #  AOS - DutyCycle, InterInstr, SpMass, Instrcount
    case3_bools = [False, False, False, False, False, False, False, False, False, True, False, False, False, False, False, False, False, True, False, False, False, False, False, False, False, True, False, False] #  AOS - InterInstr, SpMass, Instrcount
    #case3_bools = [False, False, False, False, False, True, False, False, False, True, False, False, False, False, False, False, False, True, False, False, False, False, False, False, False, True, False, False] #  AOS - Instrorb, InterInstr, SpMass, Instrcount
    #case3_bools = [True, False, False, False, True, False, False, False, True, False, False, False, False, False, False, False, True, False, False, False, True, False, False, False, True, False, False, False] #  Int Pen - DutyCycle, InstrOrb, InterInstr, SpMass, Instrsyn, Instrcount
    
    cases_dict['case1'] = case1_bools
    cases_dict['case2'] = case2_bools
    cases_dict['case3'] = case3_bools
    
    #line_colours = ['#000000','#E69F00'] # black, yellow
    line_colours = ['#000000','#E69F00','#56B4E9'] # black, yellow, blue 
    
    #casenames = ['Eps. MOEA','All heurs']
    casenames = ['Eps. MOEA','All heurs','Promising heurs']
    #casenames = ['Eps. MOEA','Bias Init - InstrCount']
    
    #alpha_values = [0.4,0.4] # change based on number of cases/visibility
    alpha_values = [0.4,0.4,0.4] # change based on number of cases/visibility
    
    line_styles = ['-','-','-']
    #line_styles = ['-','--',':']
    
    #case_hatches = ['o','/','x']
    case_hatches = [None, None, None]
else:
    case3_bools = [False, True, False, False, False, True, False, False, False, True, False, False, False, False, False, False, False, True, False, False, False, False, False, False, False, False, False, False] #  AOS - DutyCycle, InstrOrb, InterInstr, SpMass
    #case3_bools = [True, False, False, False, True, False, False, False, True, False, False, False, False, False, False, False, True, False, False, False, True, False, False, False, False, False, False, False] #  Int Pen - DutyCycle, InstrOrb, InterInstr, SpMass, Instrsyn
    
    cases_dict['case1'] = case1_bools
    cases_dict['case2'] = case2_bools
    cases_dict['case3'] = case3_bools
    
    line_colours = ['#000000','#E69F00','#56B4E9'] # black, yellow, blue
    #line_colours = ['#000000','#E69F00'] # black, yellow
    
    #casenames = ['Eps. MOEA','AOS - Heur']
    casenames = ['Eps. MOEA','All heurs','Promising heurs']
    
    line_styles = ['-','-','-']
    #line_styles = ['-','--',':']
    
    alpha_values = [0.4,0.4,0.4] # change based on number of cases/visibility
    #alpha_values = [0.5,0.5] # change based on number of cases/visibility
    
    #case_hatches = ['o','/','x']
    case_hatches = [None, None, None]

nfe_cdf_array, hv_dict_med_cases, hv_dict_1q_cases, hv_dict_3q_cases, Uvals_test, nfe_array_1 = hypervolume_computation_all_cases(cases_dict, credit_assignment, assigning_problem, num_runs, line_colours, alpha_values, casenames)

## Mann Whitney U values
print("For optimization objectives")
print(Uvals_test)

plotting_all_cases(nfe_cdf_array, hv_dict_med_cases, hv_dict_1q_cases, hv_dict_3q_cases, nfe_array_1, line_colours, alpha_values, line_styles, case_hatches, casenames)
#print(nfe_cdf_array)

plot_upto_nfe(500, 2500, hv_dict_med_cases, hv_dict_1q_cases, hv_dict_3q_cases, nfe_array_1, line_colours, alpha_values, line_styles, case_hatches, casenames)

## Checking for max NFE improvement between eps MOEA and promising heurs
hv_med_emoea = hv_dict_med_cases['case0']
hv_med_allheur = hv_dict_med_cases['case1']
hv_med_promheur = hv_dict_med_cases['case2']
hv_vals = np.linspace(0.0, 1.0, 101)
nfe_diffs_promheur = np.zeros((len(hv_vals)))
nfe_diffs_allheur = np.zeros((len(hv_vals)))
nfe_idxs_emoea = np.zeros((len(hv_vals)))
nfe_idxs_allheur = np.zeros((len(hv_vals)))
nfe_idxs_promheur = np.zeros((len(hv_vals)))
for i in range(len(hv_vals)):
    nfe_idx_emoea = find_closest_index(hv_vals[i], hv_med_emoea)
    nfe_idx_allheur = find_closest_index(hv_vals[i], hv_med_allheur)
    nfe_idx_promheur = find_closest_index(hv_vals[i], hv_med_promheur)
    nfe_diffs_promheur[i] = nfe_array[nfe_idx_emoea] - nfe_array[nfe_idx_promheur] 
    nfe_diffs_allheur[i] = nfe_array[nfe_idx_emoea] - nfe_array[nfe_idx_allheur] 
    nfe_idxs_emoea[i] = nfe_idx_emoea
    nfe_idxs_allheur[i] = nfe_idx_allheur
    nfe_idxs_promheur[i] = nfe_idx_promheur
    
max_nfe_diff_promheur = np.amax(nfe_diffs_promheur)
max_nfe_diff_allheur = np.amax(nfe_diffs_allheur)
    
