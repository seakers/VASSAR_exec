%%
clear
close
clc

%%
x = linspace(1,12,12);
lambda = 4;
f = zeros(1,size(x,2));
for i = 1:size(x,2)
    f(i) = lambda^x(i)*exp(-lambda)/factorial(x(i));
end
plot(x,f,'*b')
xlabel('k')
ylabel('Pr(X=k)')
title('Poisson distribution')