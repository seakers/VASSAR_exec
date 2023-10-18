%% Impact indices box plots
clear
close all
clc

%% Extract and store impact indices for each problem and each heuristic form
%prob_types = ['assign','partition','stiffness','artery'];
%heur_forms = ['softconstraint','operator','biasedsample'];

prob = 'artery';
heur_form = 'operator';

I_heur = read_data(prob, heur_form);

p_pos_heur = zeros(size(I_heur,2),1);
mean_heur = zeros(size(I_heur,2),1);
for i = 1:size(I_heur,2)
    p_pos_heur(i,1) = sum(I_heur(:,i) > 0)/size(I_heur,1);
    mean_heur(i,1) = mean(I_heur(:,i));
end

switch prob
    case 'assign'
        heurs = {'DC','IO','IF','PE','SM','SYN','IC'};
    case 'partition'
        heurs = {'DC','IO','IF','PE','SM','SYN'};
    otherwise
        heurs = {'PC','NP','OR','IS'};
end
heur_ticks = linspace(1, size(heurs,2), size(heurs,2));

%eoss_heurs = ['DC','IO','IF','PE','SM','SYN','IC'];
%assign_bias_heur = 'IC';
%metamat_bias_heur = 'OR';
%metamat_heurs = ['PC','NP','OR','IS'];

% Plotting
plot_boxplot(I_heur, p_pos_heur, mean_heur, heur_ticks, heurs, heur_form)

%% Plot biased sampling cases combined
probs = ["assign","stiffness","artery"];
heur_form = 'biasedsample';

heurs = {'Assignment - IC','Stiffness - OR','Artery - OR'};
heur_ticks = linspace(1, size(heurs,2), size(heurs,2));

I_heurs_combined = [];
for i = 1:size(probs,2)
    I_heur = read_data(probs(i), heur_form);
    I_heurs_combined = [I_heurs_combined, I_heur];
end

p_pos_heur = zeros(size(I_heurs_combined,2),1);
mean_heur = zeros(size(I_heurs_combined,2),1);
for i = 1:size(I_heurs_combined,2)
    p_pos_heur(i,1) = sum(I_heurs_combined(:,i) > 0)/size(I_heurs_combined,1);
    mean_heur(i,1) = mean(I_heurs_combined(:,i));
end

plot_boxplot(I_heurs_combined, p_pos_heur, mean_heur, heur_ticks, heurs, heur_form);

%% Plot truss biased sampling case

prob = "stiffness";
heur_form = 'biasedsample';

heurs = 'Stiffness - OR';
%heur_ticks = linspace(1, size(heurs,2), size(heurs,2));
heur_ticks = 1;

I_heur = read_data(prob, heur_form);

p_pos_heur = sum(I_heur > 0)/size(I_heur,1);
mean_heur = mean(I_heur);

plot_boxplot(I_heur, p_pos_heur, mean_heur, heur_ticks, heurs, heur_form);

%% Functions
function [I_heurs_array] = read_data(prob_type, heur_form_type)
    filepath = 'C:\\SEAK Lab\\SEAK Lab Github\\VASSAR\\VASSAR_exec_heur\\results\\Impact Indices\\';
    switch prob_type
        case 'assign'
            filename_prob = 'Assignment';
        case 'partition'
            filename_prob = 'Partitioning';
        case 'stiffness'
            filename_prob = 'Truss';
        case 'artery'
            filename_prob = 'Artery';
        otherwise
            filename_prob = '';
    end
    switch heur_form_type
        case 'softconstraint'
            filename_heur = 'Soft Constraint';
        case 'operator'
            filename_heur = 'Operator';
        case 'biasedsample'
            filename_heur = 'Biased Sampling';
        otherwise
            filename_heur = '';
    end
    space = ' ';
    full_filepath = [filepath,'Impact Indices - ',filename_prob,space,filename_heur,'.xlsx'];
    data_table = readtable(full_filepath);
    
    I_heurs_array = table2array(data_table(:,2:end));
end

function [] = plot_boxplot(I_heurs_arr, p_pos_heurs, mean_heurs, heur_ticks, heur_strs, heur_type)

    heur_ticks_cell = num2cell(heur_ticks);
    x_ticks = [0, heur_ticks, heur_ticks(end)+1];
    
    figure
    plot(x_ticks,zeros(size(x_ticks,2),1),'k--')
    hold on
    boxplot(I_heurs_arr, 'Labels',heur_ticks_cell)
    %hold on
    %scatter(heur_ticks, p_pos_heurs, 44, 'filled', 'pentagram', 'MarkerFaceColor', 'k')
    hold on
    scatter(heur_ticks, mean_heurs, [], 'filled', 'o', 'MarkerFaceColor', 'b')
    xlabel('Heuristic')
    ylabel('Impact Index')
    %legend('$mean(I(h))$','Interpreter','Latex','Location','Best')
    hold off
    %ylim([-1,1])
    xticks(heur_ticks)
    xticklabels(heur_strs)
    ax = gca;
    ax.FontSize = 13;
end