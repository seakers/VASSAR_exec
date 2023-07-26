%% 2-objective DTLZ1
clear
close all
clc

%% Define problem parameters
k = 5;
n_obj = 2;
n_var = n_obj + k -1;

% Variable sample points
n_samples_x = 11;
X_samples = linspace(0,1,n_samples_x);
X_des = combvec(X_samples, X_samples, X_samples, X_samples, X_samples, X_samples);

%% Compute objectives and heuristics
f_X = zeros(size(X_des, 2), n_obj);
h_X = zeros(size(X_des, 2), 1);

for i = 1:size(X_des, 2)
    X_current = X_des(:, i);
    
    % Compute objectives
    g_xm_sum = 0;
    X_k = X_current(end-k+1:end, 1);
    for j = 1:k
        g_xm_sum = g_xm_sum + ((X_k(j) - 0.5)^2 - cos(20*pi*(X_current(j) - 0.5)));
    end
    g_xm = 100*(k + g_xm_sum);
    f1_x = (1/2)*X_current(1)*(1 + g_xm);
    f2_x = (1/2)*(1 - X_current(1))*(1 + g_xm);
    f_X(i, :) = [f1_x, f2_x];
    
    % Compute heuristic
    h_X(i,1) = 1/(0.5 + sum(X_current));
end

%% Plot objectives and constraints
figure
scatter(f_X(:,1), f_X(:,2), [], h_X, 'filled')
colorbar
colormap jet
