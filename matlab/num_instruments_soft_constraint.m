%% Formulation of number of instruments heuristic - soft constraint form
clear
close all
clc

%% Function testing
x = linspace(0, 60);
model = 'logistic'; % linear or logistic

switch model
    case 'logistic'
        % Test 1 - logistic function (harsh penalization)
        L = 1;
        x0 = 30;
        k = 0.4;
        f_x = L./(1 + exp(-k*(x - x0))); 
        
    case 'linear'
        % Test 2 - Linear function (1 at 60, 0 at 15)
        f_x = zeros(1, size(x,2));
        for i = 1:size(x,2)
            if x(i) < 15
                f_x(i) = 0;
            else
                f_x(i) = 1 + (x(i) - 60)/45;
            end
        end
end

figure
plot(x, f_x)
xlabel('Number of Instruments')
ylabel('Violation Score')
