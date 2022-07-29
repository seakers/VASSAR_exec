%% Check heuristic satisfaction after application of heuristic operators
clear
close all
clc

%% Extract and analyze data for requisite problem
assigning_problem = false; % true -> assigning problem, false -> partitioning problem
move_mode = true; % true -> operators move instruments, false -> operators remove instruments (ONLY FOR ASSIGNING OPERATORS)

%filepath = "C:\\SEAK Lab\\SEAK Lab Github\\VASSAR\\VASSAR_exec_heur\\results\\";
filepath = "C:\\Users\\rosha\\Documents\\SEAK Lab Github\\VASSAR\\VASSAR_exec_heur\\results\\";

filename = "operator_heuristic_satisfaction";
if assigning_problem
    if move_mode
        filename = strcat(filename,"_assigning_move_mod.csv");
    else
        filename = strcat(filename,"_assigning_remove_mod.csv");
    end
else
    filename = strcat(filename,"_partitioning.csv");
end
filepath = strcat(filepath,filename);

format = '%s%s%f%f%s%f%f%s%f%f%s%f%f%s%f%f%s%f%f'; 
% [Full_design_initial, Full_design_instrdc, instrdc_old, instrdc_new, Full_design_instrorb, instrorb_old, instrorb_new, 
% Full_design_interinstr, interinstr_old, interinstr_new, Full_design_packeff, packeff_old, packeff_new, 
% Full_design_spmass, spmass_old, spmass_new, Full_design_instrsyn, instrsyn_old, instrsyn_new]

data_table = readtable(filepath,'Format',format,'HeaderLines',1);

instrdc_old = table2array(data_table(:,3));
instrdc_new = table2array(data_table(:,4));
instrorb_old = table2array(data_table(:,6));
instrorb_new = table2array(data_table(:,7));
interinstr_old = table2array(data_table(:,9));
interinstr_new = table2array(data_table(:,10));
packeff_old = table2array(data_table(:,12));
packeff_new = table2array(data_table(:,13));
spmass_old = table2array(data_table(:,15));
spmass_new = table2array(data_table(:,16));
instrsyn_old = table2array(data_table(:,18));
instrsyn_new = table2array(data_table(:,19));

% Peform Welch's t-test test to compare old and new datasets for different heuristics 
[~, p_instrdc] = ttest2(instrdc_old, instrdc_new, 'VarType', 'unequal', 'Tail', 'right');
[~, p_instrorb] = ttest2(instrorb_old, instrorb_new, 'VarType', 'unequal', 'Tail', 'right');
[~, p_interinstr] = ttest2(interinstr_old, interinstr_new, 'VarType', 'unequal', 'Tail', 'right');
[~, p_packeff] = ttest2(packeff_old, packeff_new, 'VarType', 'unequal', 'Tail', 'right');
[~, p_spmass] = ttest2(spmass_old, spmass_new, 'VarType', 'unequal', 'Tail', 'right');
[~, p_instrsyn] = ttest2(instrsyn_old, instrsyn_new, 'VarType', 'unequal', 'Tail', 'right');

% Plot boxplots
labels_before = repmat({'Before Operator'},size(instrdc_old,1),1);
labels_after = repmat({'After Operator'},size(instrdc_new,1),1);
figure
subplot(3,2,1)
boxplot([instrdc_old; instrdc_new],[labels_before; labels_after])
title('Duty Cycle')
subplot(3,2,2)
boxplot([instrorb_old; instrorb_new],[labels_before; labels_after])
title('Instrument Orbit Relations')
subplot(3,2,3)
boxplot([interinstr_old; interinstr_new],[labels_before; labels_after])
title('Interference')
subplot(3,2,4)
boxplot([packeff_old; packeff_new],[labels_before; labels_after])
title('Packing Efficiency')
subplot(3,2,5)
boxplot([spmass_old; spmass_new],[labels_before; labels_after])
title('Spacecraft Mass')
subplot(3,2,6)
boxplot([instrsyn_old; instrsyn_new],[labels_before; labels_after])
title('Synergy')
sgtitle('Heuristic Improvement by Operators','FontSize',10) 

%% Compute I1 using heuristic operators
I1_instrdc = mean(instrdc_new - instrdc_old);
I1_instrorb = mean(instrorb_new - instrorb_old);
I1_interinstr = mean(interinstr_new - interinstr_old);
I1_packeff = mean(packeff_new - packeff_old);
I1_spmass = mean(spmass_new - spmass_old);
I1_instrsyn = mean(instrsyn_new - instrsyn_old);

