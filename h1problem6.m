clear;
clc;
close all;

p_exact = 23*pi/192;

N_values = round(logspace(2, 7, 30));

p_estimates = zeros(size(N_values));

for i = 1:length(N_values)

    N = N_values(i);

    X = rand(N,1);
    Y = rand(N,1);
    Z = rand(N,1);

    success = (X.^2 + Y.^2 < Z) & (Z.^2 > X.*Y);

    p_estimates(i) = mean(success);

end

figure;

semilogx(N_values, p_estimates, 'o-');
hold on;

yline(p_exact, '--');

xlabel('Sample size N');
ylabel('Estimated probability');
title('Monte Carlo Estimate vs. Sample Size');
legend('Monte Carlo estimate', 'Analytic value', ...
    'Location', 'best');

grid on;