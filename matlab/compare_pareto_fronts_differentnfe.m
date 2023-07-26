%% Compare Pareto Fronts between cases
clear
close all
clc

%% Parameters
assign_prob = false;
random_init = false;

n_runs = 30;

include_case3 = true;

% Case 1 - Epsilon MOEA
case1_instrdc_bools = [false, false, false, false];
case1_instrorb_bools = [false, false, false, false];
case1_interinstr_bools = [false, false, false, false];
case1_packeff_bools = [false, false, false, false];
case1_spmass_bools = [false, false, false, false];
case1_instrsyn_bools = [false, false, false, false];
case1_heur_bools = [case1_instrdc_bools; case1_instrorb_bools; case1_interinstr_bools; case1_packeff_bools; case1_spmass_bools; case1_instrsyn_bools];

% Case 2 - AOS: All Heuristics
case2_instrdc_bools = [false, true, false, false];
case2_instrorb_bools = [false, true, false, false];
case2_interinstr_bools = [false, true, false, false];
case2_packeff_bools = [false, true, false, false];
case2_spmass_bools = [false, true, false, false];
case2_instrsyn_bools = [false, true, false, false];
case2_heur_bools = [case2_instrdc_bools; case2_instrorb_bools; case2_interinstr_bools; case2_packeff_bools; case2_spmass_bools; case2_instrsyn_bools];

% Case 3 - AOS: All Heuristics except Packing Efficiency
case3_instrdc_bools = [false, true, false, false];
case3_instrorb_bools = [false, true, false, false];
case3_interinstr_bools = [false, true, false, false];
case3_packeff_bools = [false, false, false, false];
case3_spmass_bools = [false, true, false, false];
case3_instrsyn_bools = [false, true, false, false];
case3_heur_bools = [case3_instrdc_bools; case3_instrorb_bools; case3_interinstr_bools; case3_packeff_bools; case3_spmass_bools; case3_instrsyn_bools];

%% Operation
termination_nfes = [0, 50, 100, 200, 300, 400, 500, 750, 1000, 1500];
case1_pfs = struct;
case2_pfs = struct;
if include_case3
    case3_pfs = struct;
end

for i = 1:size(termination_nfes, 2)
    nfe_thresh = termination_nfes(i);
    objs_allruns_case1 = [];
    objs_allruns_case2 = [];
    if include_case3
        objs_allruns_case3 = [];
    end
    
    for j = 1:n_runs
        [objs_runj_case1, ~, ~] = read_csv_data_tillnfe(assign_prob, case1_heur_bools, random_init, nfe_thresh, j-1);
        [objs_runj_case2, ~, ~] = read_csv_data_tillnfe(assign_prob, case2_heur_bools, random_init, nfe_thresh, j-1);
        if include_case3
            [objs_runj_case3, ~, ~] = read_csv_data_tillnfe(assign_prob, case3_heur_bools, random_init, nfe_thresh, j-1);
        end
       
        objs_allruns_case1 = vertcat(objs_allruns_case1, objs_runj_case1);
        objs_allruns_case2 = vertcat(objs_allruns_case2, objs_runj_case2);
        if include_case3
            objs_allruns_case3 = vertcat(objs_allruns_case3, objs_runj_case3);
        end
        objs_allruns_case1 = vertcat(objs_allruns_case1, objs_runj_case1);
    end
    
    current_field = strcat('nfe',num2str(nfe_thresh));
    pf_case1_nfe = compute_pareto_front(-objs_allruns_case1(:,1), objs_allruns_case1(:,2)); % science score is to be maximized 
    case1_pfs.(current_field) = pf_case1_nfe;
    
    pf_case2_nfe = compute_pareto_front(-objs_allruns_case2(:,1), objs_allruns_case2(:,2)); 
    case2_pfs.(current_field) = pf_case2_nfe;
    
    if include_case3
        pf_case3_nfe = compute_pareto_front(-objs_allruns_case3(:,1), objs_allruns_case3(:,2)); 
        case3_pfs.(current_field) = pf_case3_nfe;
    end
end

%% Plotting
for i = 1:size(termination_nfes, 2)
    nfe_thresh = termination_nfes(i);
    current_field = strcat('nfe',num2str(nfe_thresh));
    case1_pf_nfe = case1_pfs.(current_field);
    case2_pf_nfe = case2_pfs.(current_field);
    if include_case3
        case3_pf_nfe = case3_pfs.(current_field);
    end
    
    figure
    plot(-case1_pf_nfe(:,1), case1_pf_nfe(:,2), 'k*')
    hold on
    plot(-case2_pf_nfe(:,1), case2_pf_nfe(:,2), 'ro')
    
    if include_case3
        hold on
        plot(-case3_pf_nfe(:,1), case3_pf_nfe(:,2), 'bs')
        hold off
        legend('Eps. MOEA', 'AOS - All Heurs.', 'AOS - Prom. Heurs.','Location','Best')
    else
        hold off
        legend('Eps. MOEA', 'AOS - All Heurs.','Location','Best')
    end
    xlabel('Science')
    ylabel('Cost')
    title(strcat('Pareto Front Comparison - NFE = ', num2str(nfe_thresh)))
end

%% Functions
function [objs_array_req, heurs_array_req, des_array_req] = read_csv_data_tillnfe(assign_problem, heur_bools, rand_init, nfe_to_reach, run_num)
    filepath = 'C:\\SEAK Lab\\SEAK Lab Github\\VASSAR\\VASSAR_exec_heur\\results\\'; % for lab system
    %filepath = 'C:\\Users\\rosha\\Documents\\SEAK Lab Github\\VASSAR\\VASSAR_exec_heur\\results\\'; % for laptop
    methods = ["Int Pen","AOS","Bias Init","ACH"];
    heuristics = ["InstrDC","InstrOrb","InterInstr","PackEff","SpMass","InstrSyn"];
    heur_abbr = ["d","o","i","p","m","s"];
    
    filepath_random = 'Random\\';
    filename_random = 'random_';
    filename_init = '';
    filepath_credassign = 'set improvement dominance\\';
    
    eps_moea = true;
    filename_moea = 'EpsilonMOEA_emoea_';
    if (any(heur_bools(:,2)))
        filename_moea = 'AOSMOEA_emoea_';
        eps_moea = false;
    end
    constrained = '';
    filename_snippet = '';
    constr_count = 0;
    for i = 1:4
        method_heurs = '';
        method_heur_abbr = '';
        for j = 1:size(heuristics,2)
            if(heur_bools(j,i))
                method_heurs = strcat(method_heurs, heuristics(j));
                method_heur_abbr = strcat(method_heur_abbr, heur_abbr(j));
            end
        end
        if isempty(method_heurs)
            constr_count = constr_count + 1;
        else
            constrained = strcat(constrained, methods(i), " - ", method_heurs, "\\");
            filename_snippet = strcat(filename_snippet, method_heur_abbr, "con", num2str(i-1), "_");
        end
    end
    
    filepath_init = 'random initialization\\';
    if ~rand_init
        filepath_init = 'injected initialization\\';
    end
    
    filepath_moea = '';
    if (constr_count == 4)
        filepath_moea = 'Epsilon MOEA\\';
        filepath_credassign = '';
    end
    
    filepath_random = strcat(filepath_moea,constrained);
    if ~eps_moea
        filename_random = strcat(filename_moea,num2str(run_num),filename_snippet);
    else
        filename_random = strcat(filename_moea,filename_snippet);
    end
    
    filepath_prob = 'Assigning\\';
    filename_prob = '_assigning';
    if ~assign_problem
        filepath_prob = 'Partitioning\\';
        filename_prob = '_partitioning';
    end
    
    if eps_moea
        full_filename = strcat(filename_random,num2str(run_num),filename_prob,"_fullpop.csv");
    else
        full_filename = strcat(filename_random,filename_prob,"_fullpop.csv");
    end
    
    %%%% read appropriate file 
    full_filepath = strcat(filepath,filepath_prob,filepath_random,filepath_credassign,filepath_init,full_filename);
    
    format_string = '%s';
    n_data_variables = 9;

    for j = 1:n_data_variables
        format_string = strcat(format_string,'%f');
    end
        
    data_table = readtable(full_filepath,'Format',format_string,'Delimiter',',');
    
    %%%% store retrieved data into different variables
    %%%% csv_data includes: [Science, Cost, Duty Cycle Violation, Instrument Orbit Assignment Violation, 
    %%%% Interference Violation, Packing Efficiency Violation, Spacecraft Mass Violation, 
    %%%% Synergy Violation, *Instrument Count Violation* (only for Assigning Problem)]
    pop_size =  size(data_table,1);
    csv_data = zeros(pop_size,9);
    
    designs = strings(pop_size);
    designs = data_table(:,1);
        
    csv_data = data_table(:,end-8:end);

    data_array = table2array(csv_data);
    design_array = table2array(designs);
    
    full_nfe_array = data_array(:,1);
    
    init_popsize = sum(full_nfe_array == 0);
    
    [nfe_sorted, sort_indices] = sort(full_nfe_array);
    data_array_sorted = data_array(sort_indices,:);
    design_array_sorted = design_array(sort_indices,:);
    
    for i = 1:init_popsize
        nfe_sorted(i) = nfe_sorted(i) + i;
    end
    
    %nfe_sorted(init_popsize+1:end) = nfe_sorted(init_popsize+1:end) - init_popsize;
    
    %[~,closest_nfe_index] = min(abs(nfe_sorted - nfe_to_reach));
    abs_diffs = abs(nfe_sorted - nfe_to_reach);
    closest_nfe_index = find(abs_diffs == min(abs_diffs), 1, 'last');
    %nfe_array = nfe_sorted(1:closest_nfe_index);
    
    objs_array_req = data_array_sorted(1:closest_nfe_index,2:3);
    heurs_array_req = data_array_sorted(1:closest_nfe_index,4:end);
    des_array_req = design_array_sorted(1:closest_nfe_index,:);

end

function objs_pf = compute_pareto_front(objs1, objs2)
    objs = [objs1, objs2];
    pf_bool = paretofront(objs);
    objs_pf = objs(pf_bool==1,:);
end