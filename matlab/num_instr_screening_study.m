%% Visualize screening study designs for assigning problem by number of instruments
clear
close all
clc

%% Cases to consider for GA data
random_data_bool = true;
assign_case = true;

% Case 1 - Epsilon MOEA
random_init = true;
case_instrdc_bools = [false, false, false, false];
case_instrorb_bools = [false, false, false, false];
case_interinstr_bools = [false, false, false, false];
case_packeff_bools = [false, false, false, false];
case_spmass_bools = [false, false, false, false];
case_instrsyn_bools = [false, false, false, false];
case_heur_bools = [case_instrdc_bools; case_instrorb_bools; case_interinstr_bools; case_packeff_bools; case_spmass_bools; case_instrsyn_bools];

%% Read random architectures from csv
n_pop = 300; % population size of random/epsilon MOEA runs
n_runs = 10; % number of runs 

n_des = 300; % number of designs of random & epsilon MOEA each to use for correlation study

add_ga_data = true;

science_all = struct;
cost_all = struct;

n_instr_all = struct;

for i = 1:n_runs
    [f_rand_run, ~, des_rand_run] = read_csv_data(assign_case, case_heur_bools, random_data_bool, random_init, i-1);
    f_rand_red = f_rand_run(1:n_des,:);
    des_rand_red = des_rand_run(1:n_des,:);
    if add_ga_data
        [f_ga_run, ~, des_ga_run] = read_csv_data(assign_case, case_heur_bools, ~random_data_bool, random_init, i-1);
        f_ga_red = f_ga_run(1:n_des,:);
        des_ga_red = des_ga_run(1:n_des,:);
        f_rand_red = cat(1,f_rand_red,f_ga_red);
        des_rand_red = cat(1,des_rand_red,des_ga_red);
    end
    
    n_instr_run = zeros(size(des_rand_red,1),1);
    for j = 1:size(des_rand_red,1)
        n_instr_run(j,1) = get_num_instruments(des_rand_red{j});
    end
    
    current_field = strcat('trial',num2str(i));
    
    science_all.(current_field) = f_rand_red(:,1);
    cost_all.(current_field) = f_rand_red(:,2);   
    n_instr_all.(current_field) = n_instr_run;
end

%% Combining and Plotting
n_des_total = 0;
for i = 1:n_runs
    current_field = strcat('trial',num2str(i));
    science_total = science_all.(current_field);
    n_des_total = n_des_total + size(science_total,1);
end
science_array = zeros(n_des_total,1);
cost_array = zeros(n_des_total,1);
n_instr_array = zeros(n_des_total,1);
        
index = 1;
for i = 1:n_runs
    current_field = strcat('trial',num2str(i));
    science_total = science_all.(current_field);
    cost_total = cost_all.(current_field);
    n_instr_total = n_instr_all.(current_field);
    
    n_des_run = size(science_total,1);
    
    science_array(index:index+n_des_run-1,1) = science_total;
    cost_array(index:index+n_des_run-1,1) = cost_total;
    n_instr_array(index:index+n_des_run-1,1) = n_instr_total;

    index = index + n_des_run;
end

figure 
scatter(science_array,cost_array,[],n_instr_array,'filled')
xlabel('Science Score','FontSize',16)
ylabel('Cost','FontSize',16)
colorbar
title('Number of Instruments','FontSize',16)
%% Functions

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
        csv_data = data_table(:,end-7:end);
    else
        csv_data = data_table(:,end-8:end);
    end
    
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

function n_instr = get_num_instruments(arch)
    n_instr = 0;
    for i = 1:size(arch,2)
        if arch(i) == '1'
            n_instr = n_instr + 1;
        end
    end
end