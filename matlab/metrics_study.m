%% Metrics Study 
clear
close all
clc

%% Cases to consider for GA data
random_data_bool = true;

% Case 1 - Epsilon MOEA
assign_case = false;
random_init = true;
case_instrdc_bools = [false, false, false, false];
case_instrorb_bools = [false, false, false, false];
case_interinstr_bools = [false, false, false, false];
case_packeff_bools = [false, false, false, false];
case_spmass_bools = [false, false, false, false];
case_instrsyn_bools = [false, false, false, false];
case_instrcount_bools = [false, false, false, false];

if ~assign_case
    case_heur_bools = [case_instrdc_bools; case_instrorb_bools; case_interinstr_bools; case_packeff_bools; case_spmass_bools; case_instrsyn_bools];
else
    case_heur_bools = [case_instrdc_bools; case_instrorb_bools; case_interinstr_bools; case_packeff_bools; case_spmass_bools; case_instrsyn_bools; case_instrcount_bools];
end

objs_norm = [0.425, 2.5e4];
if ~assign_case
    objs_norm = [0.4, 7250];
end

%% Compute support for ease of satisfaction study
n_des2 = 300; % number of random designs to use for ease of satisfaction study
n_runs2 = 10; % number of runs 

support_rand_instrdc_runs = zeros(n_runs2,1);
support_rand_instrorb_runs = zeros(n_runs2,1);
support_rand_interinstr_runs = zeros(n_runs2,1);
support_rand_packeff_runs = zeros(n_runs2,1);
support_rand_spmass_runs = zeros(n_runs2,1);
support_rand_instrsyn_runs = zeros(n_runs2,1);
if assign_case
    support_rand_instrcount_runs = zeros(n_runs2,1);
end

for i = 1:n_runs2
    [~, heur_rand_run2, ~] = read_csv_data(assign_case, case_heur_bools, random_data_bool, random_init, i-1);
    
    heur_rand_red2 = heur_rand_run2(1:n_des2,:);
    
    instrdc_run_rand = heur_rand_red2(:,1);
    instrorb_run_rand = heur_rand_red2(:,2);
    interinstr_run_rand = heur_rand_red2(:,3);
    packeff_run_rand = heur_rand_red2(:,4);
    spmass_run_rand = heur_rand_red2(:,5);
    instrsyn_run_rand = heur_rand_red2(:,6);
    if assign_case
        instrcount_run_rand = heur_rand_red2(:,7);
    end
    
    support_rand_instrdc_runs(i,1) = length(instrdc_run_rand(instrdc_run_rand==0))/size(instrdc_run_rand,1);
    support_rand_instrorb_runs(i,1) = length(instrorb_run_rand(instrorb_run_rand==0))/size(instrdc_run_rand,1);
    support_rand_interinstr_runs(i,1) = length(interinstr_run_rand(interinstr_run_rand==0))/size(instrdc_run_rand,1);
    support_rand_packeff_runs(i,1) = length(packeff_run_rand(packeff_run_rand==0))/size(instrdc_run_rand,1);
    support_rand_spmass_runs(i,1) = length(spmass_run_rand(spmass_run_rand==0))/size(instrdc_run_rand,1);
    support_rand_instrsyn_runs(i,1) = length(instrsyn_run_rand(instrsyn_run_rand==0))/size(instrdc_run_rand,1);
    if assign_case
        support_rand_instrcount_runs(i,1) = length(instrcount_run_rand(instrcount_run_rand==0))/size(instrcount_run_rand,1);
    end
    
end

if assign_case
    support_rand_tablestats = [mean(support_rand_instrdc_runs),std(support_rand_instrdc_runs);
                            mean(support_rand_instrorb_runs),std(support_rand_instrorb_runs);
                            mean(support_rand_interinstr_runs),std(support_rand_interinstr_runs);
                            mean(support_rand_packeff_runs),std(support_rand_packeff_runs);
                            mean(support_rand_spmass_runs),std(support_rand_spmass_runs);
                            mean(support_rand_instrsyn_runs),std(support_rand_instrsyn_runs);
                            mean(support_rand_instrcount_runs),std(support_rand_instrcount_runs)];
else
    support_rand_tablestats = [mean(support_rand_instrdc_runs),std(support_rand_instrdc_runs);
                            mean(support_rand_instrorb_runs),std(support_rand_instrorb_runs);
                            mean(support_rand_interinstr_runs),std(support_rand_interinstr_runs);
                            mean(support_rand_packeff_runs),std(support_rand_packeff_runs);
                            mean(support_rand_spmass_runs),std(support_rand_spmass_runs);
                            mean(support_rand_instrsyn_runs),std(support_rand_instrsyn_runs)];
end


%% Compute objectives, constraints and heuristics for correlation study
n_pop = 300; % population size of random/epsilon MOEA runs
n_runs = 10; % number of runs 

n_des = 300; % number of designs of random & epsilon MOEA each to use for correlation study

add_ga_data = false;

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
if assign_case
    instrcount_all = struct;
end

support_full_instrdc_runs = zeros(n_runs,1);
support_full_instrorb_runs = zeros(n_runs,1);
support_full_interinstr_runs = zeros(n_runs,1);
support_full_packeff_runs = zeros(n_runs,1);
support_full_spmass_runs = zeros(n_runs,1);
support_full_instrsyn_runs = zeros(n_runs,1);
if assign_case
    support_full_instrcount_runs = zeros(n_runs,1);
end

n_pf_runs = zeros(n_runs,1);
support_pf_runs = zeros(n_runs,1);

for i = 1:n_runs
    [f_rand_run, heur_rand_run, des_rand_run] = read_csv_data(assign_case, case_heur_bools, random_data_bool, random_init, i-1);
    f_rand_red = f_rand_run(1:n_des,:).*objs_norm;
    heur_rand_red = heur_rand_run(1:n_des,:);
    des_rand_red = des_rand_run(1:n_des,:);
    if add_ga_data
        [f_ga_run, heur_ga_run, des_ga_run] = read_csv_data(assign_case, case_heur_bools, ~random_data_bool, random_init, i-1);
        f_ga_red = f_ga_run(1:n_des,:).*objs_norm;
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
    if assign_case
        instrcount_run = heur_rand_red(:,7);
    end
    
    instrdc_all.(current_field) = instrdc_run;
    instrorb_all.(current_field) = instrorb_run;
    interinstr_all.(current_field) = interinstr_run;
    packeff_all.(current_field) = packeff_run;
    spmass_all.(current_field) = spmass_run;
    instrsyn_all.(current_field) = instrsyn_run;
    if assign_case
        instrcount_all.(current_field) = instrcount_run;
    end
    
    science_true_all.(current_field) = f_rand_red(:,1);
    cost_true_all.(current_field) = f_rand_red(:,2);
    
    
    
    %%%% Normalizing objectives and pfs wrt max and min from objectives
    science_max = max(f_rand_red(:,1));
    science_min = min(f_rand_red(:,1));
    
    cost_max = max(f_rand_red(:,2));
    cost_min = min(f_rand_red(:,2));
    
    science_max_all(i,1) = science_max;
    science_min_all(i,1) = science_min;
    
    cost_max_all(i,1) = cost_max;
    cost_min_all(i,1) = cost_min;
    
    %science_norm_run = (f_rand_red(:,1) - science_min)/(science_max - science_min);
    %cost_norm_run = (f_rand_red(:,2) - cost_min)/(cost_max - cost_min);
    
    %science_all.(current_field) = science_norm_run;
    %cost_all.(current_field) = cost_norm_run;
    
    %objs_norm_pareto_run = [(objs_pareto_true(:,1) - science_min)/(science_max - science_min), (objs_pareto_true(:,2) - cost_min)/(cost_max - cost_min)]; 
    
    %min_dist_pf_run = zeros(size(f_rand_red,1),1);
    %for k = 1:size(science_norm_run,1)
        %min_dist_pf_run(k,1) = compute_min_pf_dist([science_norm_run(k,1),cost_norm_run(k,1)],objs_norm_pareto_run);
    %end
    
    %min_dist_pf_all.(current_field) = min_dist_pf_run;
    
    support_full_instrdc_runs(i,1) = length(instrdc_run(instrdc_run==0))/size(instrdc_run,1);
    support_full_instrorb_runs(i,1) = length(instrorb_run(instrorb_run==0))/size(instrdc_run,1);
    support_full_interinstr_runs(i,1) = length(interinstr_run(interinstr_run==0))/size(instrdc_run,1);
    support_full_packeff_runs(i,1) = length(packeff_run(packeff_run==0))/size(instrdc_run,1);
    support_full_spmass_runs(i,1) = length(spmass_run(spmass_run==0))/size(instrdc_run,1);
    support_full_instrsyn_runs(i,1) = length(instrsyn_run(instrsyn_run==0))/size(instrdc_run,1);
    if assign_case
        support_full_instrcount_runs(i,1) = length(instrcount_run(instrcount_run==0))/size(instrdc_run,1);
    end

    %n_pf_runs(i,1) = size(objs_pareto,1);
    
    %%%support_pf_runs(i,1) = size(objs_pareto,1)/size(science_norm_run,1);
    %support_pf_runs(i,1) = mean(min_dist_pf_run);
    
end

%%% Normalize objectives and Pareto Front
science_max_allruns = max(science_max_all);
science_min_allruns = min(science_min_all);
cost_max_allruns = max(cost_max_all);
cost_min_allruns = min(cost_min_all);
for i = 1:n_runs
    current_field = strcat('trial',num2str(i));
    
    science_run = science_true_all.(current_field);
    cost_run = cost_true_all.(current_field);
    
    objs_pareto = compute_pareto_front(-science_run, cost_run);
    objs_pareto_true = [-objs_pareto(:,1),objs_pareto(:,2)];
    
    science_norm_run = (science_run - science_min_allruns)/(science_max_allruns - science_min_allruns);
    cost_norm_run = (cost_run - cost_min_allruns)/(cost_max_allruns - cost_min_allruns);
    
    science_all.(current_field) = science_norm_run;
    cost_all.(current_field) = cost_norm_run;
    
    objs_norm_pareto_run = [(objs_pareto_true(:,1) - science_min_allruns)/(science_max_allruns - science_min_allruns), (objs_pareto_true(:,2) - cost_min_allruns)/(cost_max_allruns - cost_min_allruns)]; 
    
    min_dist_pf_run = zeros(size(science_run,1),1);
    for k = 1:size(science_norm_run,1)
        min_dist_pf_run(k,1) = compute_min_pf_dist([science_norm_run(k,1),cost_norm_run(k,1)],objs_norm_pareto_run);
    end
    
    min_dist_pf_all.(current_field) = min_dist_pf_run;
    
    n_pf_runs(i,1) = size(objs_pareto,1);
    
    %support_pf_runs(i,1) = size(objs_pareto,1)/size(science_norm_run,1);
    support_pf_runs(i,1) = mean(min_dist_pf_run);
end

if assign_case
    support_tablestats = [mean(support_pf_runs), std(support_pf_runs);
        mean(support_full_instrdc_runs), std(support_full_instrdc_runs);
        mean(support_full_instrorb_runs), std(support_full_instrorb_runs);
        mean(support_full_interinstr_runs), std(support_full_interinstr_runs);
        mean(support_full_packeff_runs), std(support_full_packeff_runs);
        mean(support_full_spmass_runs), std(support_full_spmass_runs);
        mean(support_full_instrsyn_runs), std(support_full_instrsyn_runs);
        mean(support_full_instrcount_runs), std(support_full_instrcount_runs)];
else
    support_tablestats = [mean(support_pf_runs), std(support_pf_runs);
        mean(support_full_instrdc_runs), std(support_full_instrdc_runs);
        mean(support_full_instrorb_runs), std(support_full_instrorb_runs);
        mean(support_full_interinstr_runs), std(support_full_interinstr_runs);
        mean(support_full_packeff_runs), std(support_full_packeff_runs);
        mean(support_full_spmass_runs), std(support_full_spmass_runs);
        mean(support_full_instrsyn_runs), std(support_full_instrsyn_runs)];
end

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
if assign_case
    insrcount_array = zeros(n_des_total,1);
end

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
    if assign_case
        instrcount_total = instrcount_all.(current_field);
    end
    
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
    if assign_case
        instrcount_array(index:index+n_des_run-1,1) = instrcount_total;
    end

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

instrdc_thresh_val = prctile(instrdc_array,55);
instrorb_thresh_val = prctile(instrorb_array,50);
interinstr_thresh_val = prctile(interinstr_array,55);
packeff_thresh_val = prctile(packeff_array,96);
spmass_thresh_val = prctile(spmass_array,72);
instrsyn_thresh_val = prctile(instrsyn_array,60);
if assign_case
    instrcount_thresh_val = prctile(instrcount_array,75);
end

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

if assign_case
    figure 
    scatter(science_true_array,cost_true_array,[],instrcount_array,'filled')
    xlabel('Science Score','FontSize',16)
    ylabel('Cost','FontSize',16)
    colorbar
    title('Instrument Count Violation','FontSize',16)
end

%% Compute correlation coefficients of each heuristic with objectives 
pearson_instrdc_pfdist = zeros(n_runs,1);
pearson_instrorb_pfdist = zeros(n_runs,1);
pearson_interinstr_pfdist = zeros(n_runs,1);
pearson_packeff_pfdist = zeros(n_runs,1);
pearson_spmass_pfdist = zeros(n_runs,1);
pearson_instrsyn_pfdist = zeros(n_runs,1);
if assign_case
    pearson_instrcount_pfdist = zeros(n_runs,1);
end

pval_pearson_instrdc_pfdist = zeros(n_runs,1);
pval_pearson_instrorb_pfdist = zeros(n_runs,1);
pval_pearson_interinstr_pfdist = zeros(n_runs,1);
pval_pearson_packeff_pfdist = zeros(n_runs,1);
pval_pearson_spmass_pfdist = zeros(n_runs,1);
pval_pearson_instrsyn_pfdist = zeros(n_runs,1);
if assign_case
    pval_pearson_instrcount_pfdist = zeros(n_runs,1);
end

spearman_instrdc_pfdist = zeros(n_runs,1);
spearman_instrorb_pfdist = zeros(n_runs,1);
spearman_interinstr_pfdist = zeros(n_runs,1);
spearman_packeff_pfdist = zeros(n_runs,1);
spearman_spmass_pfdist = zeros(n_runs,1);
spearman_instrsyn_pfdist = zeros(n_runs,1);
if assign_case
    spearman_instrcount_pfdist = zeros(n_runs,1);
end

pval_spearman_instrdc_pfdist = zeros(n_runs,1);
pval_spearman_instrorb_pfdist = zeros(n_runs,1);
pval_spearman_interinstr_pfdist = zeros(n_runs,1);
pval_spearman_packeff_pfdist = zeros(n_runs,1);
pval_spearman_spmass_pfdist = zeros(n_runs,1);
pval_spearman_instrsyn_pfdist = zeros(n_runs,1);
if assign_case
    pval_spearman_instrcount_pfdist = zeros(n_runs,1);
end

for i = 1:n_runs
    current_field = strcat('trial',num2str(i));
    min_dist_pf_total = min_dist_pf_all.(current_field);
    
    instrdc_total = instrdc_all.(current_field);
    instrorb_total = instrorb_all.(current_field);
    interinstr_total = interinstr_all.(current_field);
    packeff_total = packeff_all.(current_field);
    spmass_total = spmass_all.(current_field);
    instrsyn_total = instrsyn_all.(current_field);
    if assign_case
        instrcount_total = instrcount_all.(current_field);
    end
    
    [pearson_instrdc_pfdist(i),pval_pearson_instrdc_pfdist(i)] = corr(instrdc_total,min_dist_pf_total,'Type','Pearson','Rows','complete');
    [pearson_instrorb_pfdist(i),pval_pearson_instrorb_pfdist(i)] = corr(instrorb_total,min_dist_pf_total,'Type','Pearson','Rows','complete');
    [pearson_interinstr_pfdist(i),pval_pearson_interinstr_pfdist(i)] = corr(interinstr_total,min_dist_pf_total,'Type','Pearson','Rows','complete');  
    [pearson_packeff_pfdist(i),pval_pearson_packeff_pfdist(i)] = corr(packeff_total,min_dist_pf_total,'Type','Pearson','Rows','complete');
    [pearson_spmass_pfdist(i),pval_pearson_spmass_pfdist(i)] = corr(spmass_total,min_dist_pf_total,'Type','Pearson','Rows','complete');
    [pearson_instrsyn_pfdist(i),pval_pearson_instrsyn_pfdist(i)] = corr(instrsyn_total,min_dist_pf_total,'Type','Pearson','Rows','complete');
    if assign_case
        [pearson_instrcount_pfdist(i),pval_pearson_instrcount_pfdist(i)] = corr(instrcount_total,min_dist_pf_total,'Type','Pearson','Rows','complete');
    end

    [spearman_instrdc_pfdist(i),pval_spearman_instrdc_pfdist(i)] = corr(instrdc_total,min_dist_pf_total,'Type','Spearman','Rows','complete');
    [spearman_instrorb_pfdist(i),pval_spearman_instrorb_pfdist(i)] = corr(instrorb_total,min_dist_pf_total,'Type','Spearman','Rows','complete');
    [spearman_interinstr_pfdist(i),pval_spearman_interinstr_pfdist(i)] = corr(interinstr_total,min_dist_pf_total,'Type','Spearman','Rows','complete');  
    [spearman_packeff_pfdist(i),pval_spearman_packeff_pfdist(i)] = corr(packeff_total,min_dist_pf_total,'Type','Spearman','Rows','complete');
    [spearman_spmass_pfdist(i),pval_spearman_spmass_pfdist(i)] = corr(spmass_total,min_dist_pf_total,'Type','Spearman','Rows','complete');
    [spearman_instrsyn_pfdist(i),pval_spearman_instrsyn_pfdist(i)] = corr(instrsyn_total,min_dist_pf_total,'Type','Spearman','Rows','complete');
    if assign_case
        [spearman_instrcount_pfdist(i),pval_spearman_instrcount_pfdist(i)] = corr(instrcount_total,min_dist_pf_total,'Type','Spearman','Rows','complete');
    end
    
end

if assign_case
    correlation_tablestats = [mean(pearson_instrdc_pfdist), std(pearson_instrdc_pfdist), mean(spearman_instrdc_pfdist), std(spearman_instrdc_pfdist);
        mean(pearson_instrorb_pfdist), std(pearson_instrorb_pfdist), mean(spearman_instrorb_pfdist), std(spearman_instrorb_pfdist);
        mean(pearson_interinstr_pfdist), std(pearson_interinstr_pfdist), mean(spearman_interinstr_pfdist), std(spearman_interinstr_pfdist);
        mean(pearson_packeff_pfdist), std(pearson_packeff_pfdist), mean(spearman_packeff_pfdist), std(spearman_packeff_pfdist);
        mean(pearson_spmass_pfdist), std(pearson_spmass_pfdist), mean(spearman_spmass_pfdist), std(spearman_spmass_pfdist);
        mean(pearson_instrsyn_pfdist), std(pearson_instrsyn_pfdist), mean(spearman_instrsyn_pfdist), std(spearman_instrsyn_pfdist);
        mean(pearson_instrcount_pfdist), std(pearson_instrcount_pfdist), mean(spearman_instrcount_pfdist), std(spearman_instrcount_pfdist)];
    
    correlation_pval_tablestats = [mean(pval_pearson_instrdc_pfdist), std(pval_pearson_instrdc_pfdist), mean(pval_spearman_instrdc_pfdist), std(pval_spearman_instrdc_pfdist);
        mean(pval_pearson_instrorb_pfdist), std(pval_pearson_instrorb_pfdist), mean(pval_spearman_instrorb_pfdist), std(pval_spearman_instrorb_pfdist);
        mean(pval_pearson_interinstr_pfdist), std(pval_pearson_interinstr_pfdist), mean(pval_spearman_interinstr_pfdist), std(pval_spearman_interinstr_pfdist);
        mean(pval_pearson_packeff_pfdist), std(pval_pearson_packeff_pfdist), mean(pval_spearman_packeff_pfdist), std(pval_spearman_packeff_pfdist);
        mean(pval_pearson_spmass_pfdist), std(pval_pearson_spmass_pfdist), mean(pval_spearman_spmass_pfdist), std(pval_spearman_spmass_pfdist);
        mean(pval_pearson_instrsyn_pfdist), std(pval_pearson_instrsyn_pfdist), mean(pval_spearman_instrsyn_pfdist), std(pval_spearman_instrsyn_pfdist);
        mean(pval_pearson_instrcount_pfdist), std(pval_pearson_instrcount_pfdist), mean(pval_spearman_instrcount_pfdist), std(pval_spearman_instrcount_pfdist)];
else
    correlation_tablestats = [mean(pearson_instrdc_pfdist), std(pearson_instrdc_pfdist), mean(spearman_instrdc_pfdist), std(spearman_instrdc_pfdist);
        mean(pearson_instrorb_pfdist), std(pearson_instrorb_pfdist), mean(spearman_instrorb_pfdist), std(spearman_instrorb_pfdist);
        mean(pearson_interinstr_pfdist), std(pearson_interinstr_pfdist), mean(spearman_interinstr_pfdist), std(spearman_interinstr_pfdist);
        mean(pearson_packeff_pfdist), std(pearson_packeff_pfdist), mean(spearman_packeff_pfdist), std(spearman_packeff_pfdist);
        mean(pearson_spmass_pfdist), std(pearson_spmass_pfdist), mean(spearman_spmass_pfdist), std(spearman_spmass_pfdist);
        mean(pearson_instrsyn_pfdist), std(pearson_instrsyn_pfdist), mean(spearman_instrsyn_pfdist), std(spearman_instrsyn_pfdist)];
    
    correlation_pval_tablestats = [mean(pval_pearson_instrdc_pfdist), std(pval_pearson_instrdc_pfdist), mean(pval_spearman_instrdc_pfdist), std(pval_spearman_instrdc_pfdist);
        mean(pval_pearson_instrorb_pfdist), std(pval_pearson_instrorb_pfdist), mean(pval_spearman_instrorb_pfdist), std(pval_spearman_instrorb_pfdist);
        mean(pval_pearson_interinstr_pfdist), std(pval_pearson_interinstr_pfdist), mean(pval_spearman_interinstr_pfdist), std(pval_spearman_interinstr_pfdist);
        mean(pval_pearson_packeff_pfdist), std(pval_pearson_packeff_pfdist), mean(pval_spearman_packeff_pfdist), std(pval_spearman_packeff_pfdist);
        mean(pval_pearson_spmass_pfdist), std(pval_pearson_spmass_pfdist), mean(pval_spearman_spmass_pfdist), std(pval_spearman_spmass_pfdist);
        mean(pval_pearson_instrsyn_pfdist), std(pval_pearson_instrsyn_pfdist), mean(pval_spearman_instrsyn_pfdist), std(pval_spearman_instrsyn_pfdist)];
end

isnan_pearson_packeff = isnan(pearson_packeff_pfdist);
isnan_spearman_packeff = isnan(spearman_packeff_pfdist);
mean_pearson_packeff = mean(pearson_packeff_pfdist(~isnan_pearson_packeff & ~isnan_spearman_packeff));
std_pearson_packeff = std(pearson_packeff_pfdist(~isnan_pearson_packeff & ~isnan_spearman_packeff));
mean_spearman_packeff = mean(spearman_packeff_pfdist(~isnan_pearson_packeff & ~isnan_spearman_packeff));
std_spearman_packeff = std(spearman_packeff_pfdist(~isnan_pearson_packeff & ~isnan_spearman_packeff));

%% Heuristic indices calculation
corr_exp_instrdc = 1; % correlation expectation with pfdist
corr_exp_instrorb = 1; % correlation expectation with pfdist
corr_exp_interinstr = 1; % correlation expectation with pfdist
corr_exp_packeff = 1; % correlation expectation with pfdist
corr_exp_spmass = 1; % correlation expectation with pfdist
corr_exp_instrsyn = 1; % correlation expectation with pfdist
if assign_case
    corr_exp_instrcount = 1; % correlation expectation with pfdist
end
    
corr_min = 0.4; % minimum significant correlation value

%%%% VERSION 1
% indices_instrdc_total = zeros(n_runs,1);
% indices_instrorb_total = zeros(n_runs,1);
% indices_interinstr_total = zeros(n_runs,1);
% indices_packeff_total = zeros(n_runs,1);
% indices_spmass_total = zeros(n_runs,1);
% indices_instrsyn_total = zeros(n_runs,1);
% 
% indices_instrdc_norm_total = zeros(n_runs,1);
% indices_instrorb_norm_total = zeros(n_runs,1);
% indices_interinstr_norm_total = zeros(n_runs,1);
% indices_packeff_norm_total = zeros(n_runs,1);
% indices_spmass_norm_total = zeros(n_runs,1);
% indices_instrsyn_norm_total = zeros(n_runs,1);
% 
% for i = 1:n_runs
% %     indices_instrdc_total(i) = ((find_max_value_in_array([pearson_instrdc_pfdist(i),spearman_instrdc_pfdist(i)])/support_pf_runs(i))*(corr_exp_instrdc))/(support_full_instrdc_runs(i)+1e-2);
% %     indices_instrorb_total(i) = ((find_max_value_in_array([pearson_instrorb_pfdist(i),spearman_instrorb_pfdist(i)])/support_pf_runs(i))*(corr_exp_instrorb))/(support_full_instrorb_runs(i)+1e-2);
% %     indices_interinstr_total(i) = ((find_max_value_in_array([pearson_interinstr_pfdist(i),spearman_interinstr_pfdist(i)])/support_pf_runs(i))*(corr_exp_interinstr))/(support_full_interinstr_runs(i)+1e-2);
% %     indices_packeff_total(i) = ((find_max_value_in_array([pearson_packeff_pfdist(i),spearman_packeff_pfdist(i)])/support_pf_runs(i))*(corr_exp_packeff))/(support_full_packeff_runs(i)+1e-2);
% %     indices_spmass_total(i) = ((find_max_value_in_array([pearson_spmass_pfdist(i),spearman_spmass_pfdist(i)])/support_pf_runs(i))*(corr_exp_spmass))/(support_full_spmass_runs(i)+1e-2);
% %     indices_instrsyn_total(i) = ((find_max_value_in_array([pearson_instrsyn_pfdist(i),spearman_instrsyn_pfdist(i)])/support_pf_runs(i))*(corr_exp_instrsyn))/(support_full_instrsyn_runs(i)+1e-2);
% 
%     % Not incorporating support of heuristic
%     indices_instrdc_total(i) = compute_heuristic_index_contribution(pearson_instrdc_pfdist(i), spearman_instrdc_pfdist(i), corr_min, corr_exp_instrdc, support_pf_runs(i));
%     indices_instrorb_total(i) = compute_heuristic_index_contribution(pearson_instrorb_pfdist(i), spearman_instrorb_pfdist(i), corr_min, corr_exp_instrorb, support_pf_runs(i));
%     indices_interinstr_total(i) = compute_heuristic_index_contribution(pearson_interinstr_pfdist(i), spearman_interinstr_pfdist(i), corr_min, corr_exp_interinstr, support_pf_runs(i));
%     indices_packeff_total(i) = compute_heuristic_index_contribution(pearson_packeff_pfdist(i), spearman_packeff_pfdist(i), corr_min, corr_exp_packeff, support_pf_runs(i));
%     indices_spmass_total(i) = compute_heuristic_index_contribution(pearson_spmass_pfdist(i), spearman_spmass_pfdist(i), corr_min, corr_exp_spmass, support_pf_runs(i));
%     indices_instrsyn_total(i) = compute_heuristic_index_contribution(pearson_instrsyn_pfdist(i), spearman_instrsyn_pfdist(i), corr_min, corr_exp_instrsyn, support_pf_runs(i));
%     
%     indices_sum = abs(indices_instrdc_total(i)) + abs(indices_instrorb_total(i)) + abs(indices_interinstr_total(i)) + abs(indices_packeff_total(i)) + abs(indices_spmass_total(i)) + abs(indices_instrsyn_total(i));
%     indices_instrdc_norm_total(i) = indices_instrdc_total(i)/indices_sum;
%     indices_instrorb_norm_total(i)= indices_instrorb_total(i)/indices_sum;
%     indices_interinstr_norm_total(i) = indices_interinstr_total(i)/indices_sum;
%     indices_packeff_norm_total(i) = indices_packeff_total(i)/indices_sum;
%     indices_spmass_norm_total(i) = indices_spmass_total(i)/indices_sum;
%     indices_instrsyn_norm_total(i) = indices_instrsyn_total(i)/indices_sum;
% end
% 
% indices_tablestats = [mean(indices_instrdc_total),std(indices_instrdc_total);
%                       mean(indices_instrorb_total),std(indices_instrorb_total);
%                       mean(indices_interinstr_total),std(indices_interinstr_total);
%                       mean(indices_packeff_total),std(indices_packeff_total);
%                       mean(indices_spmass_total),std(indices_spmass_total);
%                       mean(indices_instrsyn_total),std(indices_instrsyn_total)];
%                   
% indices_norm_tablestats = [mean(indices_instrdc_norm_total),std(indices_instrdc_norm_total);
%                       mean(indices_instrorb_norm_total),std(indices_instrorb_norm_total);
%                       mean(indices_interinstr_norm_total),std(indices_interinstr_norm_total);
%                       mean(indices_packeff_norm_total),std(indices_packeff_norm_total);
%                       mean(indices_spmass_norm_total),std(indices_spmass_norm_total);
%                       mean(indices_instrsyn_norm_total),std(indices_instrsyn_norm_total)];
    
%%%% VERSION 2
% Compute average correlation coefficients
corr_avg_instrdc_pfdist = (pearson_instrdc_pfdist + spearman_instrdc_pfdist)/2;
corr_avg_instrorb_pfdist = (pearson_instrorb_pfdist + spearman_instrorb_pfdist)/2;
corr_avg_interinstr_pfdist = (pearson_interinstr_pfdist + spearman_interinstr_pfdist)/2;
corr_avg_packeff_pfdist = (pearson_packeff_pfdist + spearman_packeff_pfdist)/2;
corr_avg_spmass_pfdist = (pearson_spmass_pfdist + spearman_spmass_pfdist)/2;
corr_avg_instrsyn_pfdist = (pearson_instrsyn_pfdist + spearman_instrsyn_pfdist)/2;
if assign_case
    corr_avg_instrcount_pfdist = (pearson_instrcount_pfdist + spearman_instrcount_pfdist)/2;
end

% Compute I1 for each heuristic
% I1_instrdc = compute_heuristic_I1_contribution(corr_avg_instrdc_pfdist, corr_min, corr_exp_instrdc, support_pf_runs);
% I1_instrorb = compute_heuristic_I1_contribution(corr_avg_instrorb_pfdist, corr_min, corr_exp_instrorb, support_pf_runs);
% I1_interinstr = compute_heuristic_I1_contribution(corr_avg_interinstr_pfdist, corr_min, corr_exp_interinstr, support_pf_runs);
% I1_packeff = compute_heuristic_I1_contribution(corr_avg_packeff_pfdist, corr_min, corr_exp_packeff, support_pf_runs);
% I1_spmass = compute_heuristic_I1_contribution(corr_avg_spmass_pfdist, corr_min, corr_exp_spmass, support_pf_runs);
% I1_instrsyn = compute_heuristic_I1_contribution(corr_avg_instrsyn_pfdist, corr_min, corr_exp_instrsyn, support_pf_runs);

% I1_instrdc = compute_heuristic_I1_contribution2(corr_avg_instrdc_pfdist, corr_exp_instrdc, support_pf_runs);
% I1_instrorb = compute_heuristic_I1_contribution2(corr_avg_instrorb_pfdist, corr_exp_instrorb, support_pf_runs);
% I1_interinstr = compute_heuristic_I1_contribution2(corr_avg_interinstr_pfdist, corr_exp_interinstr, support_pf_runs);
% I1_packeff = compute_heuristic_I1_contribution2(corr_avg_packeff_pfdist, corr_exp_packeff, support_pf_runs);
% I1_spmass = compute_heuristic_I1_contribution2(corr_avg_spmass_pfdist, corr_exp_spmass, support_pf_runs);
% I1_instrsyn = compute_heuristic_I1_contribution2(corr_avg_instrsyn_pfdist, corr_exp_instrsyn, support_pf_runs);

% I1_instrdc = compute_heuristic_I1_contribution_run_vec(corr_avg_instrdc_pfdist, corr_exp_instrdc, support_pf_runs);
% I1_instrorb = compute_heuristic_I1_contribution_run_vec(corr_avg_instrorb_pfdist, corr_exp_instrorb, support_pf_runs);
% I1_interinstr = compute_heuristic_I1_contribution_run_vec(corr_avg_interinstr_pfdist, corr_exp_interinstr, support_pf_runs);
% I1_packeff = compute_heuristic_I1_contribution_run_vec(corr_avg_packeff_pfdist, corr_exp_packeff, support_pf_runs);
% I1_spmass = compute_heuristic_I1_contribution_run_vec(corr_avg_spmass_pfdist, corr_exp_spmass, support_pf_runs);
% I1_instrsyn = compute_heuristic_I1_contribution_run_vec(corr_avg_instrsyn_pfdist, corr_exp_instrsyn, support_pf_runs);
% if assign_case
%     I1_instrcount = compute_heuristic_I1_contribution_run_vec(corr_avg_instrcount_pfdist, corr_exp_instrcount, support_pf_runs);
% end

I1_instrdc = compute_heuristic_I1_contribution_run_vec2(corr_avg_instrdc_pfdist, corr_exp_instrdc, support_pf_runs);
I1_instrorb = compute_heuristic_I1_contribution_run_vec2(corr_avg_instrorb_pfdist, corr_exp_instrorb, support_pf_runs);
I1_interinstr = compute_heuristic_I1_contribution_run_vec2(corr_avg_interinstr_pfdist, corr_exp_interinstr, support_pf_runs);
I1_packeff = compute_heuristic_I1_contribution_run_vec2(corr_avg_packeff_pfdist, corr_exp_packeff, support_pf_runs);
I1_spmass = compute_heuristic_I1_contribution_run_vec2(corr_avg_spmass_pfdist, corr_exp_spmass, support_pf_runs);
I1_instrsyn = compute_heuristic_I1_contribution_run_vec2(corr_avg_instrsyn_pfdist, corr_exp_instrsyn, support_pf_runs);
if assign_case
    I1_instrcount = compute_heuristic_I1_contribution_run_vec2(corr_avg_instrcount_pfdist, corr_exp_instrcount, support_pf_runs);
end

%I1_norm_instrdc = I1_instrdc/(abs(I1_instrdc) + abs(I1_instrorb) + abs(I1_interinstr) + abs(I1_packeff) + abs(I1_spmass) + abs(I1_instrsyn));
%I1_norm_instrorb = I1_instrorb/(abs(I1_instrdc) + abs(I1_instrorb) + abs(I1_interinstr) + abs(I1_packeff) + abs(I1_spmass) + abs(I1_instrsyn));
%I1_norm_interinstr = I1_interinstr/(abs(I1_instrdc) + abs(I1_instrorb) + abs(I1_interinstr) + abs(I1_packeff) + abs(I1_spmass) + abs(I1_instrsyn));
%I1_norm_packeff = I1_packeff/(abs(I1_instrdc) + abs(I1_instrorb) + abs(I1_interinstr) + abs(I1_packeff) + abs(I1_spmass) + abs(I1_instrsyn));
%I1_norm_spmass = I1_spmass/(abs(I1_instrdc) + abs(I1_instrorb) + abs(I1_interinstr) + abs(I1_packeff) + abs(I1_spmass) + abs(I1_instrsyn));
%I1_norm_instrsyn = I1_instrsyn/(abs(I1_instrdc) + abs(I1_instrorb) + abs(I1_interinstr) + abs(I1_packeff) + abs(I1_spmass) + abs(I1_instrsyn));

% Compute I2 for each heuristic
% I2_instrdc = compute_heuristic_I2_contribution(corr_avg_instrdc_pfdist, corr_min, corr_exp_instrdc, support_pf_runs, support_full_instrdc_runs);
% I2_instrorb = compute_heuristic_I2_contribution(corr_avg_instrorb_pfdist, corr_min, corr_exp_instrorb, support_pf_runs, support_full_instrorb_runs);
% I2_interinstr = compute_heuristic_I2_contribution(corr_avg_interinstr_pfdist, corr_min, corr_exp_interinstr, support_pf_runs, support_full_interinstr_runs);
% I2_packeff = compute_heuristic_I2_contribution(corr_avg_packeff_pfdist, corr_min, corr_exp_packeff, support_pf_runs, support_full_packeff_runs);
% I2_spmass = compute_heuristic_I2_contribution(corr_avg_spmass_pfdist, corr_min, corr_exp_spmass, support_pf_runs, support_full_spmass_runs);
% I2_instrsyn = compute_heuristic_I2_contribution(corr_avg_instrsyn_pfdist, corr_min, corr_exp_instrsyn, support_pf_runs, support_full_instrsyn_runs);

% I2_instrdc = compute_heuristic_I2_contribution2(corr_avg_instrdc_pfdist, corr_exp_instrdc, support_pf_runs, support_full_instrdc_runs);
% I2_instrorb = compute_heuristic_I2_contribution2(corr_avg_instrorb_pfdist, corr_exp_instrorb, support_pf_runs, support_full_instrorb_runs);
% I2_interinstr = compute_heuristic_I2_contribution2(corr_avg_interinstr_pfdist, corr_exp_interinstr, support_pf_runs, support_full_interinstr_runs);
% I2_packeff = compute_heuristic_I2_contribution2(corr_avg_packeff_pfdist, corr_exp_packeff, support_pf_runs, support_full_packeff_runs);
% I2_spmass = compute_heuristic_I2_contribution2(corr_avg_spmass_pfdist, corr_exp_spmass, support_pf_runs, support_full_spmass_runs);
% I2_instrsyn = compute_heuristic_I2_contribution2(corr_avg_instrsyn_pfdist, corr_exp_instrsyn, support_pf_runs, support_full_instrsyn_runs);

I2_instrdc = compute_heuristic_I2_contribution_run_vec(corr_avg_instrdc_pfdist, corr_exp_instrdc, support_pf_runs, support_full_instrdc_runs);
I2_instrorb = compute_heuristic_I2_contribution_run_vec(corr_avg_instrorb_pfdist, corr_exp_instrorb, support_pf_runs, support_full_instrorb_runs);
I2_interinstr = compute_heuristic_I2_contribution_run_vec(corr_avg_interinstr_pfdist, corr_exp_interinstr, support_pf_runs, support_full_interinstr_runs);
I2_packeff = compute_heuristic_I2_contribution_run_vec(corr_avg_packeff_pfdist, corr_exp_packeff, support_pf_runs, support_full_packeff_runs);
I2_spmass = compute_heuristic_I2_contribution_run_vec(corr_avg_spmass_pfdist, corr_exp_spmass, support_pf_runs, support_full_spmass_runs);
I2_instrsyn = compute_heuristic_I2_contribution_run_vec(corr_avg_instrsyn_pfdist, corr_exp_instrsyn, support_pf_runs, support_full_instrsyn_runs);
if assign_case
    I2_instrcount = compute_heuristic_I2_contribution_run_vec(corr_avg_instrcount_pfdist, corr_exp_instrcount, support_pf_runs, support_full_instrcount_runs);
end

I2_abs_sum = abs(I2_instrdc) + abs(I2_instrorb) + abs(I2_interinstr) + abs(I2_packeff) + abs(I2_spmass) + abs(I2_instrsyn);
if assign_case
    I2_abs_sum = I2_abs_sum + abs(I2_instrcount);
end

I2_norm_instrdc = I2_instrdc/I2_abs_sum;
I2_norm_instrorb = I2_instrorb/I2_abs_sum;
I2_norm_interinstr = I2_interinstr/I2_abs_sum;
I2_norm_packeff = I2_packeff/I2_abs_sum;
I2_norm_spmass = I2_spmass/I2_abs_sum;
I2_norm_instrsyn = I2_instrsyn/I2_abs_sum;
if assign_case
    I2_norm_instrcount = I2_instrcount/I2_abs_sum;
end

if assign_case
    indices_tablestats = [mean(I1_instrdc), std(I1_instrdc);
                    mean(I1_instrorb), std(I1_instrorb);
                    mean(I1_interinstr), std(I1_interinstr);
                    mean(I1_packeff), std(I1_packeff);
                    mean(I1_spmass), std(I1_spmass);
                    mean(I1_instrsyn), std(I1_instrsyn);
                    mean(I1_instrcount), std(I1_instrcount)];
else
    indices_tablestats = [mean(I1_instrdc), std(I1_instrdc);
                    mean(I1_instrorb), std(I1_instrorb);
                    mean(I1_interinstr), std(I1_interinstr);
                    mean(I1_packeff), std(I1_packeff);
                    mean(I1_spmass), std(I1_spmass);
                    mean(I1_instrsyn), std(I1_instrsyn)];
end
                
%% Determine nth percentile for positive I 
if assign_case
    I_heurs_runs = [I1_instrdc, I1_instrorb, I1_interinstr, I1_packeff, I1_spmass, I1_instrsyn, I1_instrcount];
else
    I_heurs_runs = [I1_instrdc, I1_instrorb, I1_interinstr, I1_packeff, I1_spmass, I1_instrsyn];
end
%I_heurs_runs = [I1_instrdc, I1_instrorb, I1_interinstr, I1_spmass, I1_instrsyn];
%I_heurs_runs = I1_packeff;
n_percentile_heurs = zeros(size(I_heurs_runs, 2), 1);
percentile_vals = linspace(1,100,100);

for i = 1:size(I_heurs_runs, 2)
    I_currentheur = I_heurs_runs(:, i); 
    for j = 1:size(percentile_vals, 2)
        pctile = prctile(I_currentheur, percentile_vals(j));
        if (pctile > 0)
            n_percentile_heurs(i, 1) = percentile_vals(j);
            break;
        end
        if j == size(percentile_vals, 2)
            n_percentile_heurs(i, 1) = percentile_vals(j);
        end
    end
end
               
%% Thresholding heuristics, objectives and constraints into high and low (doesn't include instrcount for assigning problem)
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
function max_array_val = find_max_value_in_array(array)
    max_val = 0;
    for i1 = 1:length(array)
        if abs(array(i1)) > max_val
            max_val = array(i1);
            max_index = i1;
        end
    end
    max_array_val = array(max_index);
end

function index_contribution = compute_heuristic_index_contribution(pearson_heur_param, spearman_heur_param, min_corr_val, idx_corr_heur_param, supp_param)
    index_contribution = log10(find_max_value_in_array([pearson_heur_param,spearman_heur_param])*idx_corr_heur_param/min_corr_val)*(1/supp_param);
end

function index_contribution = compute_heuristic_I1_contribution(corr_array_heur_param, min_corr_val, idx_corr_heur_param, supp_array_param)
    % corr_array_heur_param and supp_array_param are (n x 1) arrays where n is the number of runs
    log_arg = (idx_corr_heur_param*mean(corr_array_heur_param))/(min_corr_val);
    index_contribution = log10(max([1e-4,log_arg]))*(-1*log10(mean(supp_array_param))); % Cap the log argument to 1e-4 to avoid complex values
end

function index_contribution = compute_heuristic_I2_contribution(corr_array_heur_param, min_corr_val, idx_corr_heur_param, supp_array_param, supp_array_heur)
    % corr_array_heur_param, supp_array_heur and supp_array_param are (n x 1) arrays where n is the number of runs
    log_arg = (idx_corr_heur_param*mean(corr_array_heur_param))/(min_corr_val);
    index_contribution = log10(max([1e-4,log_arg]))*(log10(mean(supp_array_param))*log10(mean(supp_array_heur))); % Cap the log argument to 1e-4 to avoid complex values
end

function index_contribution = compute_heuristic_I1_contribution2(corr_array_heur_param, idx_corr_heur_param, supp_array_param)
    % corr_array_heur_param and supp_array_param are (n x 1) arrays where n is the number of runs
    log_arg = idx_corr_heur_param*mean(corr_array_heur_param);
    index_contribution = log_arg*(-1*log10(mean(supp_array_param))); 
end

function index_contribution = compute_heuristic_I1_contribution_run_vec(corr_array_heur_param, idx_corr_heur_param, supp_array_param)
    % corr_array_heur_param and supp_array_param are (n x 1) arrays where n is the number of runs
    is_nan_vals = isnan(corr_array_heur_param); % nan values happen when all values in array are identical
    
    corr_array_heur_param = corr_array_heur_param(~is_nan_vals);
    supp_array_param = supp_array_param(~is_nan_vals);
    
    log_arg = idx_corr_heur_param.*corr_array_heur_param;
    index_contribution = log_arg.*(-1.*log10(supp_array_param));
end

function index_contribution = compute_heuristic_I1_contribution_run_vec2(corr_array_heur_param, idx_corr_heur_param, weight_array_param)
    % corr_array_heur_param and weight_array_param are (n x 1) arrays where n is the number of runs
    is_nan_vals = isnan(corr_array_heur_param); % nan values happen when all values in array are identical
    
    corr_array_heur_param = corr_array_heur_param(~is_nan_vals);
    weight_array_param = weight_array_param(~is_nan_vals);
    
    log_arg = idx_corr_heur_param.*corr_array_heur_param;
    index_contribution = log_arg.*weight_array_param;
end

function index_contribution = compute_heuristic_I2_contribution2(corr_array_heur_param, idx_corr_heur_param, supp_array_param, supp_array_heur)
    % corr_array_heur_param, supp_array_heur and supp_array_param are (n x 1) arrays where n is the number of runs
    log_arg = idx_corr_heur_param*mean(corr_array_heur_param);
    index_contribution = log_arg*(log10(mean(supp_array_param))*log10(mean(supp_array_heur))); 
end

function index_contribution = compute_heuristic_I2_contribution_run_vec(corr_array_heur_param, idx_corr_heur_param, supp_array_param, supp_array_heur)
    % corr_array_heur_param, supp_array_heur and supp_array_param are (n x 1) arrays where n is the number of runs
    log_arg = idx_corr_heur_param.*corr_array_heur_param;
    index_contribution = log_arg.*(log10(supp_array_param).*log10(supp_array_heur)); 
end

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
    filepath = 'C:\\SEAK Lab\\SEAK Lab Github\\VASSAR\\VASSAR_exec_heur\\results\\'; % for lab system
    %filepath = 'C:\\Users\\rosha\\Documents\\SEAK Lab Github\\VASSAR\\VASSAR_exec_heur\\results\\'; % for laptop
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
            filepath_moea = 'Epsilon MOEA - Metrics\\';
        end
        
        filepath_random = strcat(filepath_moea,constrained);
        filename_random = strcat(filename_moea,filename_snippet);
       
    end
    
    filepath_prob = 'Assigning\\';
    filename_prob = 'assign_';
    if ~assign_prob
        filepath_prob = 'Partitioning\\';
        filename_prob = 'partition_';
    end
    
    if random_data_read
        full_filename = strcat(filename_random,prob_name,filename_prob,filename_init,num2str(run_num),".csv");
    else
        full_filename = strcat(filename_random,prob_name,filename_prob,filename_init,num2str(run_num),"_fullpop.csv");
    end
    
    %%%% read appropriate file 
    full_filepath = strcat(filepath,filepath_random,filepath_prob,full_filename);
    
    format_string = '%s';
    if random_data_read
        n_data_variables = 8;
    else
        n_data_variables = 9;
    end
    
    if assign_prob
        n_data_variables = n_data_variables + 1; % to account for instrument count heuristic
    end

    for j = 1:n_data_variables
        format_string = strcat(format_string,'%f');
    end
        
    data_table = readtable(full_filepath,'Format',format_string,'Delimiter',',');
    
    %%%% store retrieved data into different variables
    %%%% csv_data includes: [Pen. Obj. 1, Pen.Obj. 2, True Obj. 1, True Obj. 2, Feasibility Score,
    %%%% Connectivity Score, Stiffness Ratio Constraint, Partial Collapsibility Score, 
    %%%% Nodal Properties Score, Orientation Score]
    pop_size =  size(data_table,1);
    if random_data_read
        csv_data = zeros(pop_size,8);
    else
        csv_data = zeros(pop_size,9);
    end
    
    designs = strings(pop_size);
    designs = data_table(:,1);
    
    csv_data = data_table(:,end-(n_data_variables-1):end);
    
    data_array = table2array(csv_data);
    design_array = table2array(designs);
    
    if random_data_read
        obj_array = data_array(:,1:2);
        heur_array = data_array(:,3:end);
    else
        obj_array = data_array(:,2:3);
        heur_array = data_array(:,4:end);
    end
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