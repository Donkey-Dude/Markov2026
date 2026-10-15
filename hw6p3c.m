clear; clc;

R0 = linspace(0.5, 3, 501);
Ns = [100, 1000];

figure;
hold on;

for N = Ns
    K = N/10;
    h1 = zeros(size(R0));

    for t = 1:numel(R0)
        rho = 1;
        denominator = 1;

        for j = 1:K-1
            rho = rho * N / (R0(t)*(N-j));
            denominator = denominator + rho;
        end

        h1(t) = 1 / denominator;
    end

    plot(R0, h1, 'LineWidth', 2, ...
        'DisplayName', sprintf('Exact, N = %d', N));
end

plot(R0, 1 - 1./R0, 'r--', 'LineWidth', 2, ...
    'DisplayName', '1 - 1/R_0');

yline(0, ':', 'HandleVisibility', 'off');
xline(1, ':', 'HandleVisibility', 'off');

xlabel('R_0');
ylabel('h_1');
title('Probability of reaching K = N/10 before extinction');
legend('Location', 'southeast');
grid on;
hold off;