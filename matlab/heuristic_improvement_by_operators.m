%% Check heuristic satisfaction after application of heuristic operators
clear
close all
clc

%% Extract and analyze data for requisite problem
assigning_problem = false; % true -> assigning problem, false -> partitioning problem

filepath = "C:\\SEAK Lab\\SEAK Lab Github\\VASSAR\\VASSAR_exec_heur\\results\\";
%filepath = "C:\\Users\\rosha\\Documents\\SEAK Lab Github\\VASSAR\\VASSAR_exec_heur\\results\\";

filename = "operator_heuristic_satisfaction";
if assigning_problem
    filename = strcat(filename,"_assigning.csv");
else
    filename = strcat(filename,"_partitioning.csv");
end
filepath = strcat(filepath,filename);

if assigning_problem
    format = '%s%s%f%f%s%f%f%s%f%f%s%f%f%s%f%f%s%f%f%s%f%f'; 
else
    format = '%s%s%f%f%s%f%f%s%f%f%s%f%f%s%f%f%s%f%f'; 
end
% [Full_design_initial, Full_design_instrdc, instrdc_old, instrdc_new, Full_design_instrorb, instrorb_old, instrorb_new, 
% Full_design_interinstr, interinstr_old, interinstr_new, Full_design_packeff, packeff_old, packeff_new, 
% Full_design_spmass, spmass_old, spmass_new, Full_design_instrsyn, instrsyn_old, instrsyn_new]
% assigning problem also has Full_design_instrcount, instrcount_old,
% instrcount_new

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
if assigning_problem
    instrcount_old = table2array(data_table(:,21));
    instrcount_new = table2array(data_table(:,22));
end

% Peform Welch's t-test test to compare old and new datasets for different heuristics 
[~, p_instrdc] = ttest2(instrdc_old, instrdc_new, 'VarType', 'unequal', 'Tail', 'right');
[~, p_instrorb] = ttest2(instrorb_old, instrorb_new, 'VarType', 'unequal', 'Tail', 'right');
[~, p_interinstr] = ttest2(interinstr_old, interinstr_new, 'VarType', 'unequal', 'Tail', 'right');
[~, p_packeff] = ttest2(packeff_old, packeff_new, 'VarType', 'unequal', 'Tail', 'right');
[~, p_spmass] = ttest2(spmass_old, spmass_new, 'VarType', 'unequal', 'Tail', 'right');
[~, p_instrsyn] = ttest2(instrsyn_old, instrsyn_new, 'VarType', 'unequal', 'Tail', 'right');
if assigning_problem
    [~, p_instrcount] = ttest2(instrcount_old, instrcount_new, 'VarType', 'unequal', 'Tail', 'right');
end

% Plot boxplots
n_cols = 2;
if assigning_problem
    n_cols = 3;
end
labels_before = repmat({'Before Operator'},size(instrdc_old,1),1);
labels_after = repmat({'After Operator'},size(instrdc_new,1),1);
figure
subplot(3,n_cols,1)
boxplot([instrdc_old; instrdc_new],[labels_before; labels_after])
title('Duty Cycle')
subplot(3,n_cols,2)
boxplot([instrorb_old; instrorb_new],[labels_before; labels_after])
title('Instrument Orbit Relations')
subplot(3,n_cols,3)
boxplot([interinstr_old; interinstr_new],[labels_before; labels_after])
title('Interference')
subplot(3,n_cols,4)
boxplot([packeff_old; packeff_new],[labels_before; labels_after])
title('Packing Efficiency')
subplot(3,n_cols,5)
boxplot([spmass_old; spmass_new],[labels_before; labels_after])
title('Spacecraft Mass')
subplot(3,n_cols,6)
boxplot([instrsyn_old; instrsyn_new],[labels_before; labels_after])
title('Synergy')
if assigning_problem
    subplot(3,n_cols,7)
    boxplot([instrcount_old; instrcount_new],[labels_before; labels_after])
    title('Instrument Count')
end
sgtitle('Heuristic Improvement by Operators','FontSize',10) 

%% Compute I1 using heuristic operators
I1_instrdc = mean(instrdc_new - instrdc_old);
I1_instrorb = mean(instrorb_new - instrorb_old);
I1_interinstr = mean(interinstr_new - interinstr_old);
I1_packeff = mean(packeff_new - packeff_old);
I1_spmass = mean(spmass_new - spmass_old);
I1_instrsyn = mean(instrsyn_new - instrsyn_old);
if assigning_problem 
    I1_instrcount = mean(instrcount_new - instrcount_old);
end

%% Compute Cohen's d for heuristics
n_instrdc_old = size(instrdc_old, 1);
n_instrdc_new = size(instrdc_new, 1);
n_instrorb_old = size(instrorb_old, 1);
n_instrorb_new = size(instrorb_new, 1);
n_interinstr_old = size(interinstr_old, 1);
n_interinstr_new = size(interinstr_new, 1);
n_packeff_old = size(packeff_old, 1);
n_packeff_new = size(packeff_new, 1);
n_spmass_old = size(spmass_old, 1);
n_spmass_new = size(spmass_new, 1);
n_instrsyn_old = size(instrsyn_old, 1);
n_instrsyn_new = size(instrsyn_new, 1);
if assigning_problem
    n_instrcount_old = size(instrcount_old, 1);
    n_instrcount_new = size(instrcount_new, 1);
end

s_pooled_instrdc = sqrt(((n_instrdc_old - 1)*var(instrdc_old) + (n_instrdc_new - 1)*var(instrdc_new))/(n_instrdc_old + n_instrdc_new - 2));
s_pooled_instrorb = sqrt(((n_instrorb_old - 1)*var(instrorb_old) + (n_instrorb_new - 1)*var(instrorb_new))/(n_instrorb_old + n_instrorb_new - 2));
s_pooled_interinstr = sqrt(((n_interinstr_old - 1)*var(interinstr_old) + (n_interinstr_new - 1)*var(interinstr_new))/(n_interinstr_old + n_interinstr_new - 2));
s_pooled_packeff = sqrt(((n_packeff_old - 1)*var(packeff_old) + (n_packeff_new - 1)*var(packeff_new))/(n_packeff_old + n_packeff_new - 2));
s_pooled_spmass = sqrt(((n_spmass_old - 1)*var(spmass_old) + (n_spmass_new - 1)*var(spmass_new))/(n_spmass_old + n_spmass_new - 2));
s_pooled_instrsyn = sqrt(((n_instrsyn_old - 1)*var(instrsyn_old) + (n_instrsyn_new - 1)*var(instrsyn_new))/(n_instrsyn_old + n_instrsyn_new - 2));
if assigning_problem
    s_pooled_instrcount = sqrt(((n_instrcount_old - 1)*var(instrcount_old) + (n_instrcount_new - 1)*var(instrcount_new))/(n_instrcount_old + n_instrcount_new - 2));
end

d_instrdc = (mean(instrdc_old) - mean(instrdc_new))/s_pooled_instrdc;
d_instrorb = (mean(instrorb_old) - mean(instrorb_new))/s_pooled_instrorb;
d_interinstr = (mean(interinstr_old) - mean(interinstr_new))/s_pooled_interinstr;
d_packeff = (mean(packeff_old) - mean(packeff_new))/s_pooled_packeff;
d_spmass = (mean(spmass_old) - mean(spmass_new))/s_pooled_spmass;
d_instrsyn = (mean(instrsyn_old) - mean(instrsyn_new))/s_pooled_instrsyn;
if assigning_problem
    d_instrcount = (mean(instrcount_old) - mean(instrcount_new))/s_pooled_instrcount;
end