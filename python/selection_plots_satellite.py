# -*- coding: utf-8 -*-
"""
Plot selection histories for each heuristic for the satellite problems

@author: roshan94
"""
import numpy as np
import csv
import matplotlib.pyplot as plt

assigning_problem = False

### Read the credit and quality csv files
file_loc = 'C:\\Users\\rosha\\Documents\\SEAK Lab Github\\VASSAR\\VASSAR_exec_heur\\results\\' # for laptop
#file_loc = 'C:\\SEAK Lab\\SEAK Lab Github\\VASSAR\\VASSAR_exec_heur\\results\\' # for workstation

aos_heur_bools = [True, True, True, False, True, True] # [instrdc, instrorb, interinstr, packeff, spmass, instrsyn]

credit_assign = 1 # 0 -> offspring parent dominance, 1 -> set improvement dominance, 2 -> set contribution dominance

#run_num = 1 # run number of results to read
num_runs = 30

heurs_list = ['Instrdc','Instrorb','Interinstr','Packeff','Spmass','Instrsyn']
heur_abbrvs_list = ['d','o','i','p','m','s']

filename = 'AOSMOEA_emoea_'

if assigning_problem:
    filepath_prob = 'Assigning\\'
    filename_prob = '_assigning'
else:
    filepath_prob = 'Partitioning\\'
    filename_prob = '_partitioning'
    
heurs = ''
heurs_path = 'AOS - '
for i in range(len(aos_heur_bools)):
    if aos_heur_bools[i]:
        heurs = heurs + heur_abbrvs_list[i]
        heurs_path = heurs_path + heurs_list[i] 

if all(aos_heur_bools):
    if assigning_problem:
        operator_strings = ['RepairDutyCycleAssigning+BitFlip','RepairInstrumentOrbitAssigning+BitFlip','RepairInterferenceAssigning+BitFlip','RepairPackingEfficiencyAssigning+BitFlip','RepairMassAssigning+BitFlip','RepairSynergyAssigning+BitFlip','OnePointCrossover+BitFlip']
    else:
        operator_strings = ['RepairDutyCyclePartitioning+PartitioningMutation','RepairInstrumentOrbitPartitioning+PartitioningMutation','RepairInterferencePartitioning+PartitioningMutation','RepairPackingEfficiencyPartitioning+PartitioningMutation','RepairMassPartitioning+PartitioningMutation','RepairSynergyPartitioning+PartitioningMutation','PartitioningCrossover+PartitioningMutation']
else:
    if assigning_problem:
        operator_strings = ['RepairDutyCycleAssigning+BitFlip','RepairInstrumentOrbitAssigning+BitFlip','RepairInterferenceAssigning+BitFlip','RepairMassAssigning+BitFlip','RepairSynergyAssigning+BitFlip','OnePointCrossover+BitFlip']
    else:
        operator_strings = ['RepairDutyCyclePartitioning+PartitioningMutation','RepairInstrumentOrbitPartitioning+PartitioningMutation','RepairInterferencePartitioning+PartitioningMutation','RepairMassPartitioning+PartitioningMutation','RepairSynergyPartitioning+PartitioningMutation','PartitioningCrossover+PartitioningMutation']

filepath_cred = 'offspring parent dominance\\'
if credit_assign == 1:
    filepath_cred = 'set improvement dominance\\'
elif credit_assign == 2:
    filepath_cred = 'set contribution dominance\\'

pop_size = 300
max_func_eval = 5000

nfe_array = np.linspace(pop_size+2, max_func_eval, (int((max_func_eval-(pop_size+2))/2)+1))
nfe_array = nfe_array.astype(int)

def get_csv_rows(filepath):
    rows_dict = {}
    
    with open(filepath, newline='') as csvfile:
        file_reader = csv.reader(csvfile)
        
        ind_row = 0
        for row in file_reader:
            rows_dict[ind_row] = row
            ind_row = ind_row + 1
            
    return rows_dict

def get_selection_histories_run(filename_sel, op_strings):
    sel_rows_filename = get_csv_rows(filename_sel)
    op_sels = {}
    
    for i in range(len(op_strings)):
        op_sel_nfes = []
        for j in range(len(sel_rows_filename)):
            if (sel_rows_filename[j][1] == op_strings[i]):
                op_sel_nfes.append(sel_rows_filename[j][0])
        op_sels[op_strings[i]] = op_sel_nfes
    
    return op_sels

### Get operator selection histories for each run
sel_hists_runs = {}
for i in range(num_runs):
    sel_filename = file_loc + filepath_prob + heurs_path + '\\' + filepath_cred + 'emoea_' + str(i) + heurs + 'con1_' + filename_prob + '_hist.csv'
    ops_hist_run = get_selection_histories_run(sel_filename, operator_strings)
    sel_hists_runs[i] = ops_hist_run
        
### Get selection frequencies for each heuristic across all runs
sel_freq = {}
for i in range(len(operator_strings)):
    sel_freq_op = np.zeros((len(nfe_array)))
    for j in range(len(nfe_array)):
        num_sel_op = 0
        for k in range(num_runs):
            sel_hist_run = sel_hists_runs[k]
            sel_hist_op_run = sel_hist_run[operator_strings[i]]
            if (str(nfe_array[j]) in sel_hist_op_run):
                num_sel_op = num_sel_op + 1
        sel_freq_op[j] = num_sel_op/num_runs
    sel_freq[operator_strings[i]] = sel_freq_op
        
### Plot selection frequencies as stacked bar graph
if all(aos_heur_bools):
    heur_labels = ['InstrDC','InstrOrb','InterInstr','PackEff','SpMass','InstrSyn','X+M']
    colors = ['lime','red','cyan','yellow','blue','magenta','black']
else:
    heur_labels = ['InstrDC','InstrOrb','InterInstr','SpMass','InstrSyn','X+M']
    colors = ['lime','red','cyan','blue','magenta','black']
    
fig = plt.figure()
sel_freq_sum = sel_freq[operator_strings[0]]
plt.bar(nfe_array, sel_freq_sum)
for i in range(1, len(operator_strings)):
    sel_freq_op = sel_freq[operator_strings[i]]
    plt.bar(nfe_array, sel_freq_op, bottom=sel_freq_sum, color=colors[i])
    sel_freq_sum = np.add(sel_freq_sum, sel_freq_op)
plt.xlabel('NFE',fontsize=14)
plt.ylabel('Selection Frequency',fontsize=14)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend(heur_labels, loc='upper center', bbox_to_anchor=(0.5,1.15), ncol=4, borderaxespad=0, prop={"size":12})
plt.show()