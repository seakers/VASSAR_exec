%% Observe changes in heuristic indices as a function of termination NFE for GAs
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
if assign_case
    case_instrcount_bools = [false, false, false, false];
    case_heur_bools = [case_instrdc_bools; case_instrorb_bools; case_interinstr_bools; case_packeff_bools; case_spmass_bools; case_instrsyn_bools; case_instrcount_bools];
else
    case_heur_bools = [case_instrdc_bools; case_instrorb_bools; case_interinstr_bools; case_packeff_bools; case_spmass_bools; case_instrsyn_bools];
end

%% Generate random designs
n_des = 300; % number of random designs to use for ease of satisfaction study
n_runs = 10; % number of runs 

des_rand_allruns = struct;
objs_rand_allruns = struct;
heurs_rand_allruns = struct;

for i = 1:n_runs
    [objs_rand_run, heurs_rand_run, des_rand_run] = read_csv_data_tillnfe(assign_case, case_heur_bools, random_data_bool, random_init, 5000, i-1);
    
    current_field = strcat('trial',num2str(i));
    
    objs_rand_allruns.(current_field) = objs_rand_run(1:n_des,:);
    heurs_rand_allruns.(current_field) = heurs_rand_run(1:n_des,:);
    des_rand_allruns.(current_field) = des_rand_run(1:n_des,:);
end

%% Compute heuristic indices for different termination indices for the GA runs
termination_nfes = [0, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000];
I_instrdc = zeros(size(termination_nfes, 2), n_runs);
I_instrorb = zeros(size(termination_nfes, 2), n_runs);
I_interinstr = zeros(size(termination_nfes, 2), n_runs);
I_packeff = zeros(size(termination_nfes, 2), n_runs);
I_spmass = zeros(size(termination_nfes, 2), n_runs);
I_instrsyn = zeros(size(termination_nfes, 2), n_runs);
if assign_case
    I_instrcount = zeros(size(termination_nfes, 2), n_runs);
end

for i = 1:size(termination_nfes, 2)
    I_heurs_i = compute_heuristic_indices(assign_case, case_heur_bools, objs_rand_allruns, heurs_rand_allruns, random_data_bool, random_init, termination_nfes(i), n_runs);
    I_instrdc(i,:) = I_heurs_i(:,1);
    I_instrorb(i,:) = I_heurs_i(:,2);
    I_interinstr(i,:) = I_heurs_i(:,3);
    I_packeff(i,:) = I_heurs_i(:,4);
    I_spmass(i,:) = I_heurs_i(:,5);
    I_instrsyn(i,:) = I_heurs_i(:,6);
    if assign_case
        I_instrcount(i,:) = I_heurs_i(:,7);
    end
end

% PLotting
figure
subplot(3,2,1)
errorbar(termination_nfes, mean(I_instrdc,2)', std(I_instrdc,0,2)')
hold on
plot([0, termination_nfes], zeros(1, (length(termination_nfes) +1)))
hold off
if assign_case
    ylim([-2 2])
else
    ylim([-0.25 2])
end
xlabel("Termination NFEs for GA runs")
ylabel("Heuristic Index")
title("Duty Cycle Violation")
subplot(3,2,2)
errorbar(termination_nfes, mean(I_instrorb,2)', std(I_instrorb,0,2)')
hold on
plot([0, termination_nfes], zeros(1, (length(termination_nfes) +1)))
hold off
if assign_case
    ylim([-2 2])
else
    ylim([-0.25 2])
end
xlabel("Termination NFEs for GA runs")
ylabel("Heuristic Index")
title("Instrument Orbit Relationship")
subplot(3,2,3)
errorbar(termination_nfes, mean(I_interinstr,2)', std(I_interinstr,0,2)')
hold on
plot([0, termination_nfes], zeros(1, (length(termination_nfes) +1)))
hold off
if assign_case
    ylim([-2 2])
else
    ylim([-0.25 2])
end
xlabel("Termination NFEs for GA runs")
ylabel("Heuristic Index")
title("Instrument Interference Violation")
subplot(3,2,4)
errorbar(termination_nfes, mean(I_packeff,2)', std(I_packeff,0,2)')
hold on
plot([0, termination_nfes], zeros(1, (length(termination_nfes) +1)))
hold off
if assign_case
    ylim([-2 2])
else
    ylim([-0.25 2])
end
xlabel("Termination NFEs for GA runs")
ylabel("Heuristic Index")
title("Packing Efficiency")
subplot(3,2,5)
errorbar(termination_nfes, mean(I_spmass,2)', std(I_spmass,0,2)')
hold on
plot([0, termination_nfes], zeros(1, (length(termination_nfes) +1)))
hold off
if assign_case
    ylim([-2 2])
else
    ylim([-0.25 2])
end
xlabel("Termination NFEs for GA runs")
ylabel("Heuristic Index")
title("Spacecraft Mass")
subplot(3,2,6)
errorbar(termination_nfes, mean(I_instrsyn,2)', std(I_instrsyn,0,2)')
hold on
plot([0, termination_nfes], zeros(1, (length(termination_nfes) +1)))
hold off
if assign_case
    ylim([-2 2])
else
    ylim([-0.25 2])
end
xlabel("Termination NFEs for GA runs")
ylabel("Heuristic Index")
title("Instrument Synergy Violation")

%% Functions

function [I_heurs] = compute_heuristic_indices(prob_assign_bool, heurs_bools, f_rand_allruns, heurs_rand_allruns, read_random_dataset, rand_initialize, terminate_nfe, n_runs)
    % I_heurs = [I_instrdc; I_instrorb; I_interinstr; I_packeff; I_spmass; I_instrsyn] * num_rums

    % Extract and store GA data upto the termination NFE
    [f_allgas, heurs_allgas, ~] = obtain_combined_ga_data_allruns(prob_assign_bool, heurs_bools, ~read_random_dataset, rand_initialize, terminate_nfe, n_runs);

    % Compute heuristic indices for the given termination NFE
    I_instrdc_allruns = zeros(n_runs, 1);
    I_instrorb_allruns = zeros(n_runs, 1);
    I_interinstr_allruns = zeros(n_runs, 1);
    I_packeff_allruns = zeros(n_runs, 1);
    I_spmass_allruns = zeros(n_runs, 1);
    I_instrsyn_allruns = zeros(n_runs, 1);
    if prob_assign_bool
        I_instrcount_allruns = zeros(n_runs, 1);
    end

    for i = 1:n_runs
        current_field = strcat('trial',num2str(i));

        f_rand_trial = f_rand_allruns.(current_field);
        heurs_rand_trial = heurs_rand_allruns.(current_field);

        f_currentcase = f_allgas.(current_field);
        heurs_currentcase = heurs_allgas.(current_field);
        
        f_rand_total = cat(1,f_rand_trial,f_currentcase);
        heurs_rand_total = cat(1,heurs_rand_trial,heurs_currentcase);

        instrdc_run = heurs_rand_total(:,1);
        instrorb_run = heurs_rand_total(:,2);
        interinstr_run = heurs_rand_total(:,3);
        packeff_run = heurs_rand_total(:,4);
        spmass_run = heurs_rand_total(:,5);
        instrsyn_run = heurs_rand_total(:,6);
        if prob_assign_bool
            instrcount_run = heurs_rand_total(:,7);
        end

        objs_pareto = compute_pareto_front(-f_rand_total(:,1),f_rand_total(:,2));
        objs_pareto_corrected = [-objs_pareto(:,1),objs_pareto(:,2)];

        % Normalizing objectives and pfs wrt max and min from objectives
        science_max = max(f_rand_total(:,1));
        science_min = min(f_rand_total(:,1));

        cost_max = max(f_rand_total(:,2));
        cost_min = min(f_rand_total(:,2));

        science_norm_run = (f_rand_total(:,1) - science_min)/(science_max - science_min);
        cost_norm_run = (f_rand_total(:,2) - cost_min)/(cost_max - cost_min);

        objs_norm_pareto_run = [(objs_pareto_corrected(:,1) - science_min)/(science_max - science_min), (objs_pareto_corrected(:,2) - cost_min)/(cost_max - cost_min)];

        min_dist_pf_run = zeros(size(f_rand_total,1),1);
        for k = 1:size(science_norm_run,1)
            min_dist_pf_run(k,1) = compute_min_pf_dist([science_norm_run(k,1),cost_norm_run(k,1)],objs_norm_pareto_run);
        end

%         fracsat_instrdc_runs = length(instrdc_run(instrdc_run==0))/size(science_norm_run,1);
%         fracsat_instrorb_runs = length(instrorb_run(instrorb_run==0))/size(science_norm_run,1);
%         fracsat_interinstr_runs = length(interinstr_run(interinstr_run==0))/size(science_norm_run,1);
%         fracsat_packeff_runs = length(packeff_run(packeff_run==0))/size(science_norm_run,1);
%         fracsat_spmass_runs = length(spmass_run(spmass_run==0))/size(science_norm_run,1);
%         fracsat_instrsyn_runs = length(instrsyn_run(instrsyn_run==0))/size(science_norm_run,1);
%         if prob_assign_bool
%             fracsat_instrsyn_runs = length(instrsyn_run(instrsyn_run==0))/size(science_norm_run,1);
%         end

        fracsat_pf_runs = size(objs_pareto,1)/size(science_norm_run,1);

        % Computing Pearson's Correlation Coefficients
        pearson_instrdc_pfdist = corr(instrdc_run,min_dist_pf_run,'Type','Pearson','Rows','complete');
        pearson_instrorb_pfdist= corr(instrorb_run,min_dist_pf_run,'Type','Pearson','Rows','complete');
        pearson_interinstr_pfdist = corr(interinstr_run,min_dist_pf_run,'Type','Pearson','Rows','complete');
        pearson_packeff_pfdist = corr(packeff_run,min_dist_pf_run,'Type','Pearson','Rows','complete');
        pearson_spmass_pfdist = corr(spmass_run,min_dist_pf_run,'Type','Pearson','Rows','complete');
        pearson_instrsyn_pfdist = corr(instrsyn_run,min_dist_pf_run,'Type','Pearson','Rows','complete');
        if prob_assign_bool
            pearson_instrcount_pfdist = corr(instrcount_run,min_dist_pf_run,'Type','Pearson','Rows','complete');
        end

        % Computing Spearman's Correlation Coefficients
        spearman_instrdc_pfdist = corr(instrdc_run,min_dist_pf_run,'Type','Spearman','Rows','complete');
        spearman_instrorb_pfdist = corr(instrorb_run,min_dist_pf_run,'Type','Spearman','Rows','complete');
        spearman_interinstr_pfdist = corr(interinstr_run,min_dist_pf_run,'Type','Spearman','Rows','complete');
        spearman_packeff_pfdist = corr(packeff_run,min_dist_pf_run,'Type','Spearman','Rows','complete');
        spearman_spmass_pfdist= corr(spmass_run,min_dist_pf_run,'Type','Spearman','Rows','complete');
        spearman_instrsyn_pfdist = corr(instrsyn_run,min_dist_pf_run,'Type','Spearman','Rows','complete');
        if prob_assign_bool
            spearman_instrcount_pfdist = corr(instrcount_run,min_dist_pf_run,'Type','Spearman','Rows','complete');
        end

        % Compute average correlation coefficients
        corr_avg_instrdc_pfdist = (pearson_instrdc_pfdist + spearman_instrdc_pfdist)/2;
        corr_avg_instrorb_pfdist = (pearson_instrorb_pfdist + spearman_instrorb_pfdist)/2;
        corr_avg_interinstr_pfdist = (pearson_interinstr_pfdist + spearman_interinstr_pfdist)/2;
        corr_avg_packeff_pfdist = (pearson_packeff_pfdist + spearman_packeff_pfdist)/2;
        corr_avg_spmass_pfdist = (pearson_spmass_pfdist + spearman_spmass_pfdist)/2;
        corr_avg_instrsyn_pfdist = (pearson_instrsyn_pfdist + spearman_instrsyn_pfdist)/2;
        if prob_assign_bool
            corr_avg_instrcount_pfdist = (pearson_instrcount_pfdist + spearman_instrcount_pfdist)/2;
        end

        % Computing Heuristic Indices
        corr_exp_instrdc = 1; % correlation expectation with pfdist
        corr_exp_instrorb = 1; % correlation expectation with pfdist
        corr_exp_interinstr = 1; % correlation expectation with pfdist
        corr_exp_packeff = 1; % correlation expectation with pfdist
        corr_exp_spmass = 1; % correlation expectation with pfdist
        corr_exp_instrsyn = 1; % correlation expectation with pfdist
        if prob_assign_bool
            corr_exp_instrcount = 1; % correlation expectation with pfdist
        end

        I_instrdc_allruns(i,1) = compute_heuristic_I1_contribution_run_vec(corr_avg_instrdc_pfdist, corr_exp_instrdc, fracsat_pf_runs);
        I_instrorb_allruns(i,1) = compute_heuristic_I1_contribution_run_vec(corr_avg_instrorb_pfdist, corr_exp_instrorb, fracsat_pf_runs);
        I_interinstr_allruns(i,1) = compute_heuristic_I1_contribution_run_vec(corr_avg_interinstr_pfdist, corr_exp_interinstr, fracsat_pf_runs);
        I_packeff_allruns(i,1) = compute_heuristic_I1_contribution_run_vec(corr_avg_packeff_pfdist, corr_exp_packeff, fracsat_pf_runs);
        I_spmass_allruns(i,1) = compute_heuristic_I1_contribution_run_vec(corr_avg_spmass_pfdist, corr_exp_spmass, fracsat_pf_runs);
        I_instrsyn_allruns(i,1) = compute_heuristic_I1_contribution_run_vec(corr_avg_instrsyn_pfdist, corr_exp_instrsyn, fracsat_pf_runs);
        if prob_assign_bool
            I_instrcount_allruns(i,1) = compute_heuristic_I1_contribution_run_vec(corr_avg_instrcount_pfdist, corr_exp_instrcount, fracsat_pf_runs);
        end
    end

    if prob_assign_bool
        I_heurs = [I_instrdc_allruns, I_instrorb_allruns, I_interinstr_allruns, I_packeff_allruns, I_spmass_allruns, I_instrsyn_allruns, I_instrcount_allruns];
    else
        I_heurs = [I_instrdc_allruns, I_instrorb_allruns, I_interinstr_allruns, I_packeff_allruns, I_spmass_allruns, I_instrsyn_allruns];
    end
end

function [objs_allcases, heuristics_allcases, designs_allcases] = obtain_combined_ga_data_allruns(prob_assign, heuristic_bools, read_random_data, random_init, term_nfe, num_runs)
    
    objs_allcases = struct;
    heuristics_allcases = struct;
    designs_allcases = struct;
    
    for i = 1:num_runs
        [objs_array, heurs_array, designs_array] = read_csv_data_tillnfe(prob_assign, heuristic_bools, read_random_data, random_init, term_nfe, i-1);
       
        current_field = strcat('trial',num2str(i));
        
        objs_allcases.(current_field) = objs_array;
        heuristics_allcases.(current_field) = heurs_array;
        designs_allcases.(current_field) = designs_array;
        
    end
end

function [objs_array_req, heurs_array_req, des_array_req] = read_csv_data_tillnfe(assign_prob, heur_bools, random_data_read, rand_init, nfe_to_reach, run_num)
    prob_name = 'ClimateCentric_';
    filepath = 'C:\\SEAK Lab\\SEAK Lab Github\\VASSAR\\VASSAR_exec_heur\\results\\'; % for lab system
    %filepath = 'C:\\Users\\rosha\\Documents\\SEAK Lab Github\\VASSAR\\VASSAR_exec_heur\\results\\'; % for laptop
    methods = ["Int Pen","AOS","Bias Init","ACH"];
    heuristics = ["InstrdC","Instrorb","Interinstr","Packeff","Spmass","Instrsyn","Instrcount"];
    heur_abbr = ["d","o","i","p","m","s","c"];
   
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
        n_data_variables = n_data_variables + 1;
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
        
    if random_data_read
        csv_data = data_table(:,end-(n_data_variables-1):end);
    else
        csv_data = data_table(:,end-n_data_variables:end);
    end

    data_array = table2array(csv_data);
    design_array = table2array(designs);

    if random_data_read
        objs_array_req = data_array(:,1:2);
        heurs_array_req = data_array(:,3:end);
        des_array_req = design_array;
    else
        full_nfe_array = data_array(:,1);
        [nfe_sorted, sort_indices] = sort(full_nfe_array);
        data_array_sorted = data_array(sort_indices,:);
        design_array_sorted = design_array(sort_indices,:);

        %[~,closest_nfe_index] = min(abs(nfe_sorted - nfe_to_reach));
        abs_diffs = abs(nfe_sorted - nfe_to_reach);
        closest_nfe_index = find(abs_diffs == min(abs_diffs), 1, 'last');
        %nfe_array = nfe_sorted(1:closest_nfe_index);
       
        objs_array_req = data_array_sorted(1:closest_nfe_index,2:3);
        heurs_array_req = data_array_sorted(1:closest_nfe_index,4:end);
        des_array_req = design_array_sorted(1:closest_nfe_index,:);
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

function index_contribution = compute_heuristic_I1_contribution_run_vec(corr_array_heur_param, idx_corr_heur_param, supp_array_param)
    % corr_array_heur_param and supp_array_param are (n x 1) arrays where n is the number of runs
    log_arg = idx_corr_heur_param.*corr_array_heur_param;
    index_contribution = log_arg.*(-1.*log10(supp_array_param));
end