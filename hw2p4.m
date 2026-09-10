N = 1e5;

a = 0.9;
lambda_f = 1000;
lambda_s = 10;

T = zeros(N,1);

for i = 1:N
    B = rand < a;
    U = rand;

    if B
        T(i) = -log(1-U)/lambda_f;
    else
        T(i) = -log(1-U)/lambda_s;
    end

end

empiricalMean = mean(T);

theoreticalMean = a/lambda_f + (1-a)/lambda_s;

threshold = 0.050;

empiricalProb = mean(T > threshold);

theoreticalProb = ...
    a*exp(-lambda_f*threshold) ...
    + (1-a)*exp(-lambda_s*threshold);

fprintf('Empirical mean:   %.6f s\n', empiricalMean);
fprintf('Theoretical mean: %.6f s\n\n', theoreticalMean);

fprintf('Empirical P(T > 50 ms):   %.6f\n', empiricalProb);
fprintf('Theoretical P(T > 50 ms): %.6f\n', theoreticalProb);

binWidth = 0.001;

histogram(T, ...
    'BinWidth', binWidth, ...
    'Normalization', 'pdf');

hold on

t = linspace(0, max(T), 2000);

f = a*lambda_f*exp(-lambda_f*t) ...
    + (1-a)*lambda_s*exp(-lambda_s*t);

plot(t, f, 'LineWidth', 2)

set(gca, 'YScale', 'log')

xlabel('Dwell Time (s)')
ylabel('Probability Density')
title(sprintf('Exponential Mixture, N = %d, Bin Width = %.3f s', ...
    N, binWidth))

legend('Empirical PDF', 'Theoretical PDF')
grid on

hold off