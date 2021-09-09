%% Metrics Study 
clear
close all
clc

%% Cases to consider for GA data
random_data_bool = true;
% Case 1 - Epsilon MOEA
assign_case = true;
random_init = true;
case_instrdc_bools = [false, false, false, false];
case_instrorb_bools = [false, false, false, false];
case_interinstr_bools = [false, false, false, false];
case_packeff_bools = [false, false, false, false];
case_spmass_bools = [false, false, false, false];
case_instrsyn_bools = [false, false, false, false];
case_heur_bools = [case_instrdc_bools; case_instrorb_bools; case_interinstr_bools; case_packeff_bools; case_spmass_bools; case_instrsyn_bools];

%% Compute objectives, constraints and heuristics
n_pop = 300; % population size of random/epsilon MOEA runs
n_runs = 10; % number of runs 

n_des = 300; % number of designs of random & epsilon MOEA each to use for metrics computation

add_ga_data = true;

min_dist_pf_all = struct;

science_all = struct;
cost_all = struct;

science_true_all = struct;
cost_true_all = struct;

science_max_all = zeros(n_runs,1);
science_min_all = zeros(n_runs,1);
cost_max_all = zeros(n_runs,1);
cost_min_all = zeros(n_runs,1);

instrdc_all = struct;
instrorb_all = struct;
interinstr_all = struct;
packeff_all = struct;
spmass_all = struct;
instrsyn_all = struct;

support_full_instrdc_runs = zeros(n_runs,1);
support_full_instrorb_runs = zeros(n_runs,1);
support_full_interinstr_runs = zeros(n_runs,1);
support_full_packeff_runs = zeros(n_runs,1);
support_full_spmass_runs = zeros(n_runs,1);
support_full_instrsyn_runs = zeros(n_runs,1);

n_pf_runs = zeros(n_runs,1);
support_pf_runs = zeros(n_runs,1);

for i = 1:n_runs
    [f_rand_run, heur_rand_run, des_rand_run] = read_csv_data(assign_case, case_heur_bools, random_data_bool, random_init, i-1);
    f_rand_red = f_rand_run(1:n_des,:);
    heur_rand_red = heur_rand_run(1:n_des,:);
    des_rand_red = des_rand_run(1:n_des,:);
    if add_ga_data
        [f_ga_run, heur_ga_run, des_ga_run] = read_csv_data(assign_case, case_heur_bools, ~random_data_bool, random_init, i-1);
        f_ga_red = f_ga_run(1:n_des,:);
        heur_ga_red = heur_ga_run(1:n_des,:);
        des_ga_red = des_ga_run(1:n_des,:);
        
        f_rand_red = cat(1,f_rand_red,f_ga_red);
        heur_rand_red = cat(1,heur_rand_red,heur_ga_red);
        des_rand_red = cat(1,des_rand_red,des_ga_red);
    end
    
    current_field = strcat('trial',num2str(i));
    
    instrdc_run = heur_rand_red(:,1);
    instrorb_run = heur_rand_red(:,2);
    interinstr_run = heur_rand_red(:,3);
    packeff_run = heur_rand_red(:,4);
    spmass_run = heur_rand_red(:,5);
    instrsyn_run = heur_rand_red(:,6);
    
    instrdc_all.(current_field) = instrdc_run;
    instrorb_all.(current_field) = instrorb_run;
    interinstr_all.(current_field) = interinstr_run;
    packeff_all.(current_field) = packeff_run;
    spmass_all.(current_field) = spmass_run;
    instrsyn_all.(current_field) = instrsyn_run;
    
    science_true_all.(current_field) = f_rand_red(:,1);
    cost_true_all.(current_field) = f_rand_red(:,2);
    
    objs_pareto = compute_pareto_front(-f_rand_red(:,1),f_rand_red(:,2));
    objs_pareto_true = [-objs_pareto(:,1),objs_pareto(:,2)];
    
    %%%% Normalizing objectives and pfs wrt max and min from objectives
    science_max = max(f_rand_red(:,1));
    science_min = min(f_rand_red(:,1));
    
    cost_max = max(f_rand_red(:,2));
    cost_min = min(f_rand_red(:,2));
    
    science_max_all(i,1) = science_max;
    science_min_all(i,1) = science_min;
    
    cost_max_all(i,1) = cost_max;
    cost_min_all(i,1) = cost_min;
    
    science_norm_run = (f_rand_red(:,1) - science_min)/(science_max - science_min);
    cost_norm_run = (f_rand_red(:,2) - cost_min)/(cost_max - cost_min);
    
    science_all.(current_field) = science_norm_run;
    cost_all.(current_field) = cost_norm_run;
    
    objs_norm_pareto_run = [(objs_pareto_true(:,1) - science_min)/(science_max - science_min), (objs_pareto_true(:,2) - cost_min)/(cost_max - cost_min)]; 
    
    min_dist_pf_run = zeros(size(f_rand_red,1),1);
    for k = 1:size(science_norm_run,1)
        min_dist_pf_run(k,1) = compute_min_pf_dist([science_norm_run(k,1),cost_norm_run(k,1)],objs_norm_pareto_run);
    end
    
    min_dist_pf_all.(current_field) = min_dist_pf_run;
    
    support_full_instrdc_runs(i,1) = length(instrdc_run(instrdc_run==0))/size(science_norm_run,1);
    support_full_instrorb_runs(i,1) = length(instrorb_run(instrorb_run==0))/size(science_norm_run,1);
    support_full_interinstr_runs(i,1) = length(interinstr_run(interinstr_run==0))/size(science_norm_run,1);
    support_full_packeff_runs(i,1) = length(packeff_run(packeff_run==0))/size(science_norm_run,1);
    support_full_spmass_runs(i,1) = length(spmass_run(spmass_run==0))/size(science_norm_run,1);
    support_full_instrsyn_runs(i,1) = length(instrsyn_run(instrsyn_run==0))/size(science_norm_run,1);

    n_pf_runs(i,1) = size(objs_pareto,1);
    
    support_pf_runs(i,1) = size(objs_pareto,1)/size(science_norm_run,1);
    
end
support_tablestats = [mean(support_pf_runs), std(support_pf_runs);
    mean(support_full_instrdc_runs), std(support_full_instrdc_runs);
    mean(support_full_instrorb_runs), std(support_full_instrorb_runs);
    mean(support_full_interinstr_runs), std(support_full_interinstr_runs);
    mean(support_full_packeff_runs), std(support_full_packeff_runs);
    mean(support_full_spmass_runs), std(support_full_spmass_runs);
    mean(support_full_instrsyn_runs), std(support_full_instrsyn_runs)];

%% Heuristic violation plots

n_des_total = 0;
for i = 1:n_runs
    current_field = strcat('trial',num2str(i));
    science_total = science_all.(current_field);
    n_des_total = n_des_total + size(science_total,1);
end
science_array = zeros(n_des_total,1);
cost_array = zeros(n_des_total,1);

science_true_array = zeros(n_des_total,1);
cost_true_array = zeros(n_des_total,1);

insrtdc_array = zeros(n_des_total,1);
instrorb_array = zeros(n_des_total,1);
interinstr_array = zeros(n_des_total,1);
packeff_array = zeros(n_des_total,1);
spmass_array = zeros(n_des_total,1);
instrsyn_array = zeros(n_des_total,1);

min_dist_pf_array = zeros(n_des_total,1);
        
index = 1;
for i = 1:n_runs
    current_field = strcat('trial',num2str(i));
    science_total = science_all.(current_field);
    cost_total = cost_all.(current_field);
    
    science_true_total = science_true_all.(current_field);
    cost_true_total = cost_true_all.(current_field);
    
    instrdc_total = instrdc_all.(current_field);
    instrorb_total = instrorb_all.(current_field);
    interinstr_total = interinstr_all.(current_field);
    packeff_total = packeff_all.(current_field);
    spmass_total = spmass_all.(current_field);
    instrsyn_total = instrsyn_all.(current_field);
    
    min_dist_pf_total = min_dist_pf_all.(current_field);
    
    n_des_run = size(science_total,1);
    
    science_array(index:index+n_des_run-1,1) = science_total;
    cost_array(index:index+n_des_run-1,1) = cost_total;
    
    science_true_array(index:index+n_des_run-1,1) = science_true_total;
    cost_true_array(index:index+n_des_run-1,1) = cost_true_total;
    
    instrdc_array(index:index+n_des_run-1,1) = instrdc_total;
    instrorb_array(index:index+n_des_run-1,1) = instrorb_total;
    interinstr_array(index:index+n_des_run-1,1) = interinstr_total;
    packeff_array(index:index+n_des_run-1,1) = packeff_total;
    spmass_array(index:index+n_des_run-1,1) = spmass_total;
    instrsyn_array(index:index+n_des_run-1,1) = instrsyn_total;

    min_dist_pf_array(index:index+n_des_run-1,1) = min_dist_pf_total;

    index = index + n_des_run;
end

% Set threshold limit
science_thresh_val = prctile(science_array,75);
cost_thresh_val = prctile(cost_array,25);

science_thresh_val_true = zeros(n_runs,1);
cost_thresh_val_true = zeros(n_runs,1);
for i = 1:n_runs
    science_thresh_val_true(i,1) = science_thresh_val*(science_max_all(i,1) - science_min_all(i,1)) + science_min_all(i,1);
    cost_thresh_val_true(i,1) = cost_thresh_val*(cost_max_all(i,1) - cost_min_all(i,1)) + cost_min_all(i,1);
end

instrdc_thresh_val = prctile(instrdc_array,25);
instrorb_thresh_val = prctile(instrorb_array,50);
interinstr_thresh_val = prctile(interinstr_array,50);
packeff_thresh_val = prctile(packeff_array,50);
spmass_thresh_val = prctile(spmass_array,50);
instrsyn_thresh_val = prctile(instrsyn_array,25);

mindist_pf_thresh_val = prctile(min_dist_pf_array,50);

score_tabelstats = [mean(instrdc_array), std(instrdc_array);
    mean(instrorb_array), std(instrorb_array);
    mean(interinstr_array), std(interinstr_array);
    mean(packeff_array), std(packeff_array);
    mean(spmass_array), std(spmass_array);
    mean(instrsyn_array), std(instrsyn_array)];

% Plot Heuristics Violations
figure 
scatter(science_true_array,cost_true_array,[],instrdc_array,'filled')
xlabel('Science Score','FontSize',16)
ylabel('Cost','FontSize',16)
colorbar
title('Instrument Duty Cycle Violation','FontSize',16)

figure 
scatter(science_true_array,cost_true_array,[],instrorb_array,'filled')
xlabel('Science Score','FontSize',16)
ylabel('Cost','FontSize',16)
colorbar
title('Instrument Orbit Relationship Violation','FontSize',16)

figure 
scatter(science_true_array,cost_true_array,[],interinstr_array,'filled')
xlabel('Science Score','FontSize',16)
ylabel('Cost','FontSize',16)
colorbar
title('Interfering Instruments Violation','FontSize',16)

figure 
scatter(science_true_array,cost_true_array,[],packeff_array,'filled')
xlabel('Science Score','FontSize',16)
ylabel('Cost','FontSize',16)
colorbar
title('Packing Efficiency Violation','FontSize',16)

figure 
scatter(science_true_array,cost_true_array,[],spmass_array,'filled')
xlabel('Science Score','FontSize',16)
ylabel('Cost','FontSize',16)
colorbar
title('Spacecraft Mass Violation','FontSize',16)

figure 
scatter(science_true_array,cost_true_array,[],instrsyn_array,'filled')
xlabel('Science Score','FontSize',16)
ylabel('Cost','FontSize',16)
colorbar
title('Instrument Synergy Violation','FontSize',16)

%% Compute correlation coefficients of each heuristic with objectives 
pearson_instrdc_pfdist = zeros(n_runs,1);
pearson_instrorb_pfdist = zeros(n_runs,1);
pearson_interinstr_pfdist = zeros(n_runs,1);
pearson_packeff_pfdist = zeros(n_runs,1);
pearson_spmass_pfdist = zeros(n_runs,1);
pearson_instrsyn_pfdist = zeros(n_runs,1);

spearman_instrdc_pfdist = zeros(n_runs,1);
spearman_instrorb_pfdist = zeros(n_runs,1);
spearman_interinstr_pfdist = zeros(n_runs,1);
spearman_packeff_pfdist = zeros(n_runs,1);
spearman_spmass_pfdist = zeros(n_runs,1);
spearman_instrsyn_pfdist = zeros(n_runs,1);

for i = 1:n_runs
    current_field = strcat('trial',num2str(i));
    min_dist_pf_total = min_dist_pf_all.(current_field);
    
    instrdc_total = instrdc_all.(current_field);
    instrorb_total = instrorb_all.(current_field);
    interinstr_total = interinstr_all.(current_field);
    packeff_total = packeff_all.(current_field);
    spmass_total = spmass_all.(current_field);
    instrsyn_total = instrsyn_all.(current_field);
    
    [pearson_instrdc_pfdist(i),~] = corr(instrdc_total,min_dist_pf_total,'Type','Pearson','Rows','complete');
    [pearson_instrorb_pfdist(i),~] = corr(instrorb_total,min_dist_pf_total,'Type','Pearson','Rows','complete');
    [pearson_interinstr_pfdist(i),~] = corr(interinstr_total,min_dist_pf_total,'Type','Pearson','Rows','complete');  
    [pearson_packeff_pfdist(i),~] = corr(packeff_total,min_dist_pf_total,'Type','Pearson','Rows','complete');
    [pearson_spmass_pfdist(i),~] = corr(spmass_total,min_dist_pf_total,'Type','Pearson','Rows','complete');
    [pearson_instrsyn_pfdist(i),~] = corr(instrsyn_total,min_dist_pf_total,'Type','Pearson','Rows','complete');

    [spearman_instrdc_pfdist(i),~] = corr(instrdc_total,min_dist_pf_total,'Type','Spearman','Rows','complete');
    [spearman_instrorb_pfdist(i),~] = corr(instrorb_total,min_dist_pf_total,'Type','Spearman','Rows','complete');
    [spearman_interinstr_pfdist(i),~] = corr(interinstr_total,min_dist_pf_total,'Type','Spearman','Rows','complete');  
    [spearman_packeff_pfdist(i),~] = corr(packeff_total,min_dist_pf_total,'Type','Spearman','Rows','complete');
    [spearman_spmass_pfdist(i),~] = corr(spmass_total,min_dist_pf_total,'Type','Spearman','Rows','complete');
    [spearman_instrsyn_pfdist(i),~] = corr(instrsyn_total,min_dist_pf_total,'Type','Spearman','Rows','complete');
    
end

correlation_tablestats = [mean(pearson_instrdc_pfdist), std(pearson_instrdc_pfdist), mean(spearman_instrdc_pfdist), std(spearman_instrdc_pfdist);
    mean(pearson_instrorb_pfdist), std(pearson_instrorb_pfdist), mean(spearman_instrorb_pfdist), std(spearman_instrorb_pfdist);
    mean(pearson_interinstr_pfdist), std(pearson_interinstr_pfdist), mean(spearman_interinstr_pfdist), std(spearman_interinstr_pfdist);
    mean(pearson_packeff_pfdist), std(pearson_packeff_pfdist), mean(spearman_packeff_pfdist), std(spearman_packeff_pfdist);
    mean(pearson_spmass_pfdist), std(pearson_spmass_pfdist), mean(spearman_spmass_pfdist), std(spearman_spmass_pfdist);
    mean(pearson_instrsyn_pfdist), std(pearson_instrsyn_pfdist), mean(spearman_instrsyn_pfdist), std(spearman_instrsyn_pfdist)];

%% Thresholding heuristics, objectives and constraints into high and low
instrdc_all_thresh = struct;
instrorb_all_thresh = struct;
interinstr_all_thresh = struct;
packeff_all_thresh = struct;
spmass_all_thresh = struct;
instrsyn_all_thresh = struct;

science_all_thresh = struct;
cost_all_thresh = struct;

min_dist_pf_all_thresh = struct;

for i = 1:n_runs
    current_field = strcat('trial',num2str(i));
    
    instrdc_total = instrdc_all.(current_field);
    instrorb_total = instrorb_all.(current_field);
    interinstr_total = interinstr_all.(current_field);
    packeff_total = packeff_all.(current_field);
    spmass_total = spmass_all.(current_field);
    instrsyn_total = instrsyn_all.(current_field);
    
    science_total = science_all.(current_field);
    cost_total = cost_all.(current_field);
    
    min_dist_pf_total = min_dist_pf_all.(current_field);
    
    instrdc_run_thresh = double(instrdc_total >= instrdc_thresh_val);
    instrorb_run_thresh = double(instrorb_total >= instrorb_thresh_val);
    interinstr_run_thresh = double(interinstr_total >= interinstr_thresh_val);
    packeff_run_thresh = double(packeff_total >= packeff_thresh_val);
    spmass_run_thresh = double(spmass_total >= spmass_thresh_val);
    instrsyn_run_thresh = double(instrsyn_total >= instrsyn_thresh_val);
    
    science_run_thresh = double(science_total >= science_thresh_val);
    cost_run_thresh = double(cost_total >= cost_thresh_val);
    
    min_dist_pf_thresh = double(min_dist_pf_total >= mindist_pf_thresh_val);
    
    instrdc_all_thresh.(current_field) = instrdc_run_thresh;
    instrorb_all_thresh.(current_field) = instrorb_run_thresh;
    interinstr_all_thresh.(current_field) = interinstr_run_thresh;
    packeff_all_thresh.(current_field) = packeff_run_thresh;
    spmass_all_thresh.(current_field) = spmass_run_thresh;
    instrsyn_all_thresh.(current_field) = instrsyn_run_thresh;
    
    science_all_thresh.(current_field) = science_run_thresh;
    cost_all_thresh.(current_field) = cost_run_thresh;
    
    min_dist_pf_all_thresh.(current_field) = min_dist_pf_thresh;
end

%% Computing interestingness measures for different association rules
% intrness = [support_a, support_b, support_ab, confidence_a2b,
% confidence_b2a, lift]
% heur_intrness_arrays_run = [low_heur_intrness; high_heur_intrness]{size = 4 x 6 x 2} one 3D array for each run  

instrdc_intrness_truepf_arrays = struct;
instrorb_intrness_truepf_arrays = struct;
interinstr_intrness_truepf_arrays = struct;
packeff_intrness_truepf_arrays = struct;
spmass_intrness_truepf_arrays = struct;
instrsyn_intrness_truepf_arrays = struct;

instrdc_intrness_truepf = struct;
instrorb_intrness_truepf = struct;
interinstr_intrness_truepf = struct;
packeff_intrness_truepf = struct;
spmass_intrness_truepf = struct;
instrsyn_intrness_truepf = struct;

for i = 1:n_runs
    
    current_field = strcat('trial',num2str(i));
    science_thresh = science_all_thresh.(current_field);
    cost_thresh = cost_all_thresh.(current_field);
    
    instrdc_thresh = instrdc_all_thresh.(current_field);
    instrorb_thresh = instrorb_all_thresh.(current_field);
    interinstr_thresh = interinstr_all_thresh.(current_field);
    packeff_thresh = packeff_all_thresh.(current_field);
    spmass_thresh = spmass_all_thresh.(current_field);
    instrsyn_thresh = instrsyn_all_thresh.(current_field);
    
    min_dist_pf_thresh = min_dist_pf_all_thresh.(current_field);
    
    instrdc_intrness_truepf_run = compute_heur_intrness_truepf(instrdc_thresh, min_dist_pf_thresh);
    instrdc_intrness_truepf.(current_field) = instrdc_intrness_truepf_run;
    instrorb_intrness_truepf_run = compute_heur_intrness_truepf(instrorb_thresh, min_dist_pf_thresh);
    instrorb_intrness_truepf.(current_field) = instrorb_intrness_truepf_run;
    interinstr_intrness_truepf_run = compute_heur_intrness_truepf(interinstr_thresh, min_dist_pf_thresh);
    interinstr_intrness_truepf.(current_field) = interinstr_intrness_truepf_run;
    packeff_intrness_truepf_run = compute_heur_intrness_truepf(packeff_thresh, min_dist_pf_thresh);
    packeff_intrness_truepf.(current_field) = packeff_intrness_truepf_run;
    spmass_intrness_truepf_run = compute_heur_intrness_truepf(spmass_thresh, min_dist_pf_thresh);
    spmass_intrness_truepf.(current_field) = spmass_intrness_truepf_run;
    instrsyn_intrness_truepf_run = compute_heur_intrness_truepf(instrsyn_thresh, min_dist_pf_thresh);
    instrsyn_intrness_truepf.(current_field) = instrsyn_intrness_truepf_run;
   
    instrdc_intrness_run = compute_heur_intrness_array(instrdc_thresh, min_dist_pf_thresh, science_thresh, cost_thresh);
    instrdc_intrness_truepf_arrays.(current_field) = instrdc_intrness_run;
    instrorb_intrness_run = compute_heur_intrness_array(instrorb_thresh, min_dist_pf_thresh, science_thresh, cost_thresh);
    instrorb_intrness_truepf_arrays.(current_field) = instrorb_intrness_run;
    interinstr_intrness_run = compute_heur_intrness_array(interinstr_thresh, min_dist_pf_thresh, science_thresh, cost_thresh);
    interinstr_intrness_truepf_arrays.(current_field) = interinstr_intrness_run;
    packeff_intrness_run = compute_heur_intrness_array(packeff_thresh, min_dist_pf_thresh, science_thresh, cost_thresh);
    packeff_intrness_truepf_arrays.(current_field) = packeff_intrness_run;
    spmass_intrness_run = compute_heur_intrness_array(spmass_thresh, min_dist_pf_thresh, science_thresh, cost_thresh);
    spmass_intrness_truepf_arrays.(current_field) = spmass_intrness_run;
    instrsyn_intrness_run = compute_heur_intrness_array(instrsyn_thresh, min_dist_pf_thresh, science_thresh, cost_thresh);
    instrsyn_intrness_truepf_arrays.(current_field) = instrsyn_intrness_run;
    
end

%% Obtain interestingness measure arrays for each heuristic in turns
low_heur_close2truepf_lowscience_runs = zeros(n_runs,6);
low_heur_close2truepf_highscience_runs = zeros(n_runs,6);
low_heur_close2truepf_lowcost_runs = zeros(n_runs,6);
low_heur_close2truepf_highcost_runs = zeros(n_runs,6);

high_heur_close2truepf_lowscience_runs = zeros(n_runs,6);
high_heur_close2truepf_highscience_runs = zeros(n_runs,6);
high_heur_close2truepf_lowcost_runs = zeros(n_runs,6);
high_heur_close2truepf_highcost_runs = zeros(n_runs,6);

low_heur_close2truepf = zeros(n_runs,6);
high_heur_close2truepf = zeros(n_runs,6);

for i = 1:n_runs
    current_field = strcat('trial',num2str(i));
    %heur_intrness_run = instrdc_intrness_truepf_arrays.(current_field); % size: 4 x 6 x 2
    %heur_intrness_run = instrorb_intrness_truepf_arrays.(current_field);
    %heur_intrness_run = interinstr_intrness_truepf_arrays.(current_field);
    %heur_intrness_run = packeff_intrness_truepf_arrays.(current_field); 
    %heur_intrness_run = spmass_intrness_truepf_arrays.(current_field);
    heur_intrness_run = instrsyn_intrness_truepf_arrays.(current_field); 
    
    low_heur_close2truepf_lowscience_runs(i,:) = heur_intrness_run(1,:,1);
    low_heur_close2truepf_highscience_runs(i,:) = heur_intrness_run(2,:,1); 
    low_heur_close2truepf_lowcost_runs(i,:) = heur_intrness_run(3,:,1);
    low_heur_close2truepf_highcost_runs(i,:) = heur_intrness_run(4,:,1);
    
    high_heur_close2truepf_lowscience_runs(i,:) = heur_intrness_run(1,:,2);
    high_heur_close2truepf_highscience_runs(i,:) = heur_intrness_run(2,:,2);
    high_heur_close2truepf_lowcost_runs(i,:) = heur_intrness_run(3,:,2);
    high_heur_close2truepf_highcost_runs(i,:) = heur_intrness_run(4,:,2);
    
    %heur_intrness_truepf_run = instrdc_intrness_truepf.(current_field); % size: 2 x 6
    %heur_intrness_truepf_run = instrorb_intrness_truepf.(current_field);
    %heur_intrness_truepf_run = interinstr_intrness_truepf.(current_field);
    %heur_intrness_truepf_run = packeff_intrness_truepf.(current_field);
    %heur_intrness_truepf_run = spmass_intrness_truepf.(current_field);
    heur_intrness_truepf_run = instrsyn_intrness_truepf.(current_field);
    
    low_heur_close2truepf(i,:) = heur_intrness_truepf_run(1,:);
    high_heur_close2truepf(i,:) = heur_intrness_truepf_run(2,:);
    
end

low_heur_tablestats = [mean(low_heur_close2truepf_lowscience_runs(:,1)), std(low_heur_close2truepf_lowscience_runs(:,1)), mean(low_heur_close2truepf_lowscience_runs(:,2)), std(low_heur_close2truepf_lowscience_runs(:,2)), mean(low_heur_close2truepf_lowscience_runs(:,3)), std(low_heur_close2truepf_lowscience_runs(:,3)), mean(low_heur_close2truepf_lowscience_runs(:,4)), std(low_heur_close2truepf_lowscience_runs(:,4)), mean(low_heur_close2truepf_lowscience_runs(:,5)), std(low_heur_close2truepf_lowscience_runs(:,5)), mean(low_heur_close2truepf_lowscience_runs(:,6)), std(low_heur_close2truepf_lowscience_runs(:,6));
    mean(low_heur_close2truepf_highscience_runs(:,1)), std(low_heur_close2truepf_highscience_runs(:,1)), mean(low_heur_close2truepf_highscience_runs(:,2)), std(low_heur_close2truepf_highscience_runs(:,2)), mean(low_heur_close2truepf_highscience_runs(:,3)), std(low_heur_close2truepf_highscience_runs(:,3)), mean(low_heur_close2truepf_highscience_runs(:,4)), std(low_heur_close2truepf_highscience_runs(:,4)), mean(low_heur_close2truepf_highscience_runs(:,5)), std(low_heur_close2truepf_highscience_runs(:,5)), mean(low_heur_close2truepf_highscience_runs(:,6)), std(low_heur_close2truepf_highscience_runs(:,6));
    mean(low_heur_close2truepf_lowcost_runs(:,1)), std(low_heur_close2truepf_lowcost_runs(:,1)), mean(low_heur_close2truepf_lowcost_runs(:,2)), std(low_heur_close2truepf_lowcost_runs(:,2)), mean(low_heur_close2truepf_lowcost_runs(:,3)), std(low_heur_close2truepf_lowcost_runs(:,3)), mean(low_heur_close2truepf_lowcost_runs(:,4)), std(low_heur_close2truepf_lowcost_runs(:,4)), mean(low_heur_close2truepf_lowcost_runs(:,5)), std(low_heur_close2truepf_lowcost_runs(:,5)), mean(low_heur_close2truepf_lowcost_runs(:,6)), std(low_heur_close2truepf_lowcost_runs(:,6));
    mean(low_heur_close2truepf_highcost_runs(:,1)), std(low_heur_close2truepf_highcost_runs(:,1)), mean(low_heur_close2truepf_highcost_runs(:,2)), std(low_heur_close2truepf_highcost_runs(:,2)), mean(low_heur_close2truepf_highcost_runs(:,3)), std(low_heur_close2truepf_highcost_runs(:,3)), mean(low_heur_close2truepf_highcost_runs(:,4)), std(low_heur_close2truepf_highcost_runs(:,4)), mean(low_heur_close2truepf_highcost_runs(:,5)), std(low_heur_close2truepf_highcost_runs(:,5)), mean(low_heur_close2truepf_highcost_runs(:,6)), std(low_heur_close2truepf_highcost_runs(:,6))];

high_heur_tablestats = [mean(high_heur_close2truepf_lowscience_runs(:,1)), std(high_heur_close2truepf_lowscience_runs(:,1)), mean(high_heur_close2truepf_lowscience_runs(:,2)), std(high_heur_close2truepf_lowscience_runs(:,2)), mean(high_heur_close2truepf_lowscience_runs(:,3)), std(high_heur_close2truepf_lowscience_runs(:,3)), mean(high_heur_close2truepf_lowscience_runs(:,4)), std(high_heur_close2truepf_lowscience_runs(:,4)), mean(high_heur_close2truepf_lowscience_runs(:,5)), std(high_heur_close2truepf_lowscience_runs(:,5)), mean(high_heur_close2truepf_lowscience_runs(:,6)), std(high_heur_close2truepf_lowscience_runs(:,6));
    mean(high_heur_close2truepf_highscience_runs(:,1)), std(high_heur_close2truepf_highscience_runs(:,1)), mean(high_heur_close2truepf_highscience_runs(:,2)), std(high_heur_close2truepf_highscience_runs(:,2)), mean(high_heur_close2truepf_highscience_runs(:,3)), std(high_heur_close2truepf_highscience_runs(:,3)), mean(high_heur_close2truepf_highscience_runs(:,4)), std(high_heur_close2truepf_highscience_runs(:,4)), mean(high_heur_close2truepf_highscience_runs(:,5)), std(high_heur_close2truepf_highscience_runs(:,5)), mean(high_heur_close2truepf_highscience_runs(:,6)), std(high_heur_close2truepf_highscience_runs(:,6));
    mean(high_heur_close2truepf_lowcost_runs(:,1)), std(high_heur_close2truepf_lowcost_runs(:,1)), mean(high_heur_close2truepf_lowcost_runs(:,2)), std(high_heur_close2truepf_lowcost_runs(:,2)), mean(high_heur_close2truepf_lowcost_runs(:,3)), std(high_heur_close2truepf_lowcost_runs(:,3)), mean(high_heur_close2truepf_lowcost_runs(:,4)), std(high_heur_close2truepf_lowcost_runs(:,4)), mean(high_heur_close2truepf_lowcost_runs(:,5)), std(high_heur_close2truepf_lowcost_runs(:,5)), mean(high_heur_close2truepf_lowcost_runs(:,6)), std(high_heur_close2truepf_lowcost_runs(:,6));
    mean(high_heur_close2truepf_highcost_runs(:,1)), std(high_heur_close2truepf_highcost_runs(:,1)), mean(high_heur_close2truepf_highcost_runs(:,2)), std(high_heur_close2truepf_highcost_runs(:,2)), mean(high_heur_close2truepf_highcost_runs(:,3)), std(high_heur_close2truepf_highcost_runs(:,3)), mean(high_heur_close2truepf_highcost_runs(:,4)), std(high_heur_close2truepf_highcost_runs(:,4)), mean(high_heur_close2truepf_highcost_runs(:,5)), std(high_heur_close2truepf_highcost_runs(:,5)), mean(high_heur_close2truepf_highcost_runs(:,6)), std(high_heur_close2truepf_highcost_runs(:,6))];

heur_truepf_tablestats = [mean(low_heur_close2truepf(:,1)), std(low_heur_close2truepf(:,1)), mean(low_heur_close2truepf(:,2)), std(low_heur_close2truepf(:,2)), mean(low_heur_close2truepf(:,3)), std(low_heur_close2truepf(:,3)), mean(low_heur_close2truepf(:,4)), std(low_heur_close2truepf(:,4)), mean(low_heur_close2truepf(:,5)), std(low_heur_close2truepf(:,5)), mean(low_heur_close2truepf(:,6)), std(low_heur_close2truepf(:,6));
    mean(high_heur_close2truepf(:,1)), std(high_heur_close2truepf(:,1)), mean(high_heur_close2truepf(:,2)), std(high_heur_close2truepf(:,2)), mean(high_heur_close2truepf(:,3)), std(high_heur_close2truepf(:,3)), mean(high_heur_close2truepf(:,4)), std(high_heur_close2truepf(:,4)), mean(high_heur_close2truepf(:,5)), std(high_heur_close2truepf(:,5)), mean(high_heur_close2truepf(:,6)), std(high_heur_close2truepf(:,6))];

%% Functions
function plot_regression(linmodel,X,y,Xname,Yname,run_num)
    X_norm = (X - mean(X))/std(X);
    Y_pred_norm = predict(linmodel,X_norm);
    Y_pred = (Y_pred_norm*std(y)) + mean(y);
    figure
    plot(X, y, '*b')
    hold on
    plot(X, Y_pred, 'k')
    hold off
    xlabel(Xname,'Interpreter','Latex')
    ylabel(Yname,'Interpreter','Latex')
    title(strcat('Linear Regression fit for run ',num2str(run_num)))
end

function heur_intrness_array = compute_heur_intrness_array(heur_thresh, mindist_pf_thresh, sciencethresh, costthresh)
    heur_intrness_array = zeros(4,6,2);
    
    [sup_lowheur_close2pf_lowscience, conf_lowheur_close2pf_lowscience, lift_lowheur_close2pf_lowscience] = compute_interestingness_measures(sum(~heur_thresh), sum(~mindist_pf_thresh & ~sciencethresh), sum(~heur_thresh & ~mindist_pf_thresh & ~sciencethresh), size(sciencethresh,1));
    [sup_lowheur_close2pf_highscience, conf_lowheur_close2pf_highscience, lift_lowheur_close2pf_highscience] = compute_interestingness_measures(sum(~heur_thresh), sum(~mindist_pf_thresh & sciencethresh), sum(~heur_thresh & ~mindist_pf_thresh & sciencethresh), size(sciencethresh,1));
    [sup_lowheur_close2pf_lowcost, conf_lowheur_close2pf_lowcost, lift_lowheur_close2pf_lowcost] = compute_interestingness_measures(sum(~heur_thresh), sum(~mindist_pf_thresh & ~costthresh), sum(~heur_thresh & ~mindist_pf_thresh & ~costthresh), size(sciencethresh,1));
    [sup_lowheur_close2pf_highcost, conf_lowcoll_close2pf_highcost, lift_lowheur_close2pf_highcost] = compute_interestingness_measures(sum(~heur_thresh), sum(~mindist_pf_thresh & costthresh), sum(~heur_thresh & ~mindist_pf_thresh & costthresh), size(sciencethresh,1));
    
    [sup_highheur_close2pf_lowscience, conf_highheur_close2pf_lowscience, lift_highheur_close2pf_lowscience] = compute_interestingness_measures(sum(heur_thresh), sum(~mindist_pf_thresh & ~sciencethresh), sum(heur_thresh & ~mindist_pf_thresh & ~sciencethresh), size(sciencethresh,1));
    [sup_highheur_close2pf_highscience, conf_highheur_close2pf_highscience, lift_highheur_close2pf_highscience] = compute_interestingness_measures(sum(heur_thresh), sum(~mindist_pf_thresh & sciencethresh), sum(heur_thresh & ~mindist_pf_thresh & sciencethresh), size(sciencethresh,1));
    [sup_highheur_close2pf_lowcost, conf_highheur_close2pf_lowcost, lift_highheur_close2pf_lowcost] = compute_interestingness_measures(sum(heur_thresh), sum(~mindist_pf_thresh & ~costthresh), sum(heur_thresh & ~mindist_pf_thresh & ~costthresh), size(sciencethresh,1));
    [sup_highheur_close2pf_highcost, conf_highheur_close2pf_highcost, lift_highheur_close2pf_highcost] = compute_interestingness_measures(sum(heur_thresh), sum(~mindist_pf_thresh & costthresh), sum(heur_thresh & ~mindist_pf_thresh & costthresh), size(sciencethresh,1));
    
    heur_intrness_array(:,:,1) = [sup_lowheur_close2pf_lowscience, conf_lowheur_close2pf_lowscience, lift_lowheur_close2pf_lowscience;
                            sup_lowheur_close2pf_highscience, conf_lowheur_close2pf_highscience, lift_lowheur_close2pf_highscience;
                            sup_lowheur_close2pf_lowcost, conf_lowheur_close2pf_lowcost, lift_lowheur_close2pf_lowcost;
                            sup_lowheur_close2pf_highcost, conf_lowcoll_close2pf_highcost, lift_lowheur_close2pf_highcost];
                        
     heur_intrness_array(:,:,2) = [sup_highheur_close2pf_lowscience, conf_highheur_close2pf_lowscience, lift_highheur_close2pf_lowscience;
                            sup_highheur_close2pf_highscience, conf_highheur_close2pf_highscience, lift_highheur_close2pf_highscience;
                            sup_highheur_close2pf_lowcost, conf_highheur_close2pf_lowcost, lift_highheur_close2pf_lowcost;
                            sup_highheur_close2pf_highcost, conf_highheur_close2pf_highcost, lift_highheur_close2pf_highcost];
                        
end

function [support_vals, confidence_vals, lift_val] = compute_interestingness_measures(n_a,n_b,n_ab,n_all)
    % compute support, confidence and lift vals for association rule A->B
    sup_a = compute_support_arm(n_a,n_all);
    sup_b = compute_support_arm(n_b,n_all);
    sup_ab = compute_support_arm(n_ab,n_all);
    conf_a2b = compute_confidence_arm(n_ab,n_a);
    conf_b2a = compute_confidence_arm(n_ab,n_b);
    lift_val = compute_lift_arm(n_a,n_b,n_ab,n_all);
    support_vals = [sup_a, sup_b, sup_ab];
    confidence_vals = [conf_a2b, conf_b2a];
    
end

function heur_intr_array = compute_heur_intrness_truepf(heur_thresh, mindist_pf_thresh)
    heur_intr_array = zeros(2,6);
    [sup_lowheur_close2pf, conf_lowheur_close2pf, lift_lowheur_close2pf] = compute_interestingness_measures(sum(~heur_thresh), sum(~mindist_pf_thresh), sum(~heur_thresh & ~mindist_pf_thresh), size(mindist_pf_thresh,1));
    [sup_highheur_close2pf, conf_highheur_close2pf, lift_highheur_close2pf] = compute_interestingness_measures(sum(heur_thresh), sum(~mindist_pf_thresh), sum(heur_thresh & ~mindist_pf_thresh), size(mindist_pf_thresh,1));
    
    heur_intr_array(1,:) = [sup_lowheur_close2pf, conf_lowheur_close2pf, lift_lowheur_close2pf];
    heur_intr_array(2,:) = [sup_highheur_close2pf, conf_highheur_close2pf, lift_highheur_close2pf];
end

function lift = compute_lift_arm(n_X,n_Y,n_XY,n_total) 
    lift = (n_XY/n_X)/(n_Y/n_total);
end

function support = compute_support_arm(n_X,n_total)
    support = n_X/n_total;
end

function confidence = compute_confidence_arm(n_XY,n_X)
    confidence = n_XY/n_X;
end

function [obj_array, heur_array, design_array] = read_csv_data(assign_prob, heur_bools, random_data_read, rand_init, run_num)
    prob_name = 'ClimateCentric_';
    filepath = 'C:\\SEAK Lab\\SEAK Lab Github\\VASSAR\\VASSAR_exec_heur\\results\\';
    methods = ["Int Pen","AOS","Bias Init","ACH"];
    heuristics = ["InstrDC","InstrOrb","InterInstr","PackEff","SpMass","InstrSyn"];
    heur_abbr = ["d","o","i","p","m","s"];
   
    filepath_random = 'Random\\';
    filename_random = 'random_';
    filename_init = '';
    if ~random_data_read
        filename_moea = 'EpsilonMOEA_emoea_';
        if (any(heur_bools(:,2)))
            filename_moea = 'AOSMOEA_emoea_';    
        end
        constrained = '';
        filename_snippet = '';
        constr_count = 0;
        for i = 1:4
            method_heurs = '';
            method_heur_abbr = '';
            for j = 1:size(heuristics,2)
                if(heur_bools(j,i))
                    method_heurs = strcat(method_heurs, heuristics(j), " ");
                    method_heur_abbr = strcat(method_heur_abbr, heur_abbr(j));
                end
            end
            if isempty(method_heurs)
                constr_count = constr_count + 1;
            else
                constrained = strcat(constrained, methods(i), " - ", method_heurs, "\\");
                filename_snippet = strcat(filename_snippet, method_heur_abbr, num2str(i), "_");
            end
        end
        
        filename_init = 'randinit_';
        if ~rand_init
            filename_init = 'randinjinit_';
        end
        
        filepath_moea = '';
        if (constr_count == 4)
            filepath_moea = 'Epsilon MOEA\\';
        end
        
        filepath_random = strcat(filepath_moea,constrained);
        filename_random = strcat(filename_moea,filename_snippet);
       
    end
    
    filepath_prob = 'Assigning\\';
    filename_prob = 'assign_';
    if ~assign_prob
        filepath_prob = 'Partitioning and Assigning\\';
        filename_prob = 'partandassign_';
    end
    
    full_filename = strcat(filename_random,prob_name,filename_prob,filename_init,num2str(run_num),".csv");
    
    %%%% read appropriate file 
    full_filepath = strcat(filepath,filepath_random,filepath_prob,full_filename);
    
    if assign_prob
        format_string = '%s';
        for j = 1:8
            format_string = strcat(format_string,'%f');
        end
    else
        format_string = '';
        for j = 1:(24+8) % 2*number of instruments in ClimateCentric params
            format_string = strcat(format_string,'%f');
        end
    end
        
    data_table = readtable(full_filepath,'Format',format_string,'HeaderLines',1);
    
    %%%% store retrieved data into different variables
    %%%% csv_data includes: [Pen. Obj. 1, Pen.Obj. 2, True Obj. 1, True Obj. 2, Feasibility Score,
    %%%% Connectivity Score, Stiffness Ratio Constraint, Partial Collapsibility Score, 
    %%%% Nodal Properties Score, Orientation Score]
    pop_size =  size(data_table,1);
    csv_data = zeros(pop_size,8);
    
    if assign_prob
        designs = strings(pop_size);
        designs = data_table(:,1);
    else
        designs = zeros(pop_size,24);
        designs = data_table(:,1:24);
    end
        
    csv_data = data_table(:,end-7:end);
    
    data_array = table2array(csv_data);
    design_array = table2array(designs);
    
    obj_array = data_array(:,1:2);
    heur_array = data_array(:,3:end);
end

function norm_array = normalize_array(array)
    norm_array = zeros(size(array,1),size(array,2));
    for i = 1:size(array,2)
        array_col = array(:,i);
        max_array_col = max(array_col);
        min_array_col = min(array_col);
        norm_array_col = (array_col - min_array_col)/(max_array_col - min_array_col);
        norm_array(:,i) = norm_array_col;
    end
end

function objs_pf = compute_pareto_front(obj1_run, obj2_run)
    objs = [obj1_run, obj2_run];
    pf_bool = paretofront(objs);
    objs_pf = objs(pf_bool==1,:);
end

function min_dist_pf = compute_min_pf_dist(x_vals_norm, pf_objs_norm)
    dist_vals = zeros(size(pf_objs_norm,1),1);
    for i = 1:size(pf_objs_norm,1)
        dist_vals(i) = sqrt((x_vals_norm(1) - pf_objs_norm(i,1))^2 + (x_vals_norm(2) - pf_objs_norm(i,2))^2);
    end
    min_dist_pf = min(dist_vals);
end