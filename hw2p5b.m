N = 1e4;
R = 1e3;
mu = 1;

C = ones(R,1);

for n = 3:N-1
    p = mu*C/n;
    grows = rand(R,1) < p;
    C = C + grows;
end

z = C/N;

empiricalMean = mean(z);
empiricalRatio = std(z)/empiricalMean;

theoreticalMean = 1/3;
theoreticalRatio = 1/sqrt(2);

fprintf('Empirical mean of z: %.4f\n', empiricalMean);
fprintf('Theoretical mean:    %.4f\n', theoreticalMean);
fprintf('Empirical SD/mean:   %.4f\n', empiricalRatio);
fprintf('Theoretical SD/mean: %.4f\n', theoreticalRatio);
fprintf('Smallest core: %d\n', min(C));
fprintf('Largest core:  %d\n', max(C));

binWidth = 0.05;

histogram(z, ...
    'BinWidth', binWidth, ...
    'Normalization', 'pdf');

hold on

x = linspace(0,1,1000);
h = 2*(1-x);

plot(x, h, 'LineWidth', 2)

xlabel('z = C/N')
ylabel('Probability Density')
title(sprintf('Core Growth Process, N = %d, R = %d', N, R))
legend('Empirical PDF', 'h(z) = 2(1-z)')

hold off