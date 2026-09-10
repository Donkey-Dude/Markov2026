N = 1e4;
binWidth = 0.2;

samples = zeros(N,1);
accepted = 0;
proposals = 0;
lambda = 0.2;
c = 1/(0.16*exp(1));

tic

while accepted < N
    U1 = rand;
    X = -log(1-U1)/lambda;
    U2 = rand;
    proposals = proposals + 1;
    if U2 < (5/c)*X*exp(-0.8*X)
        accepted = accepted + 1;
        samples(accepted) = X;
    end

end

totalTime = toc;

empiricalAcceptance = accepted/proposals;

theoreticalAcceptance = 1/c;

meanTime = totalTime/N;

fprintf('Number accepted: %d\n', accepted);
fprintf('Number proposed: %d\n', proposals);
fprintf('Empirical acceptance fraction: %.4f\n', empiricalAcceptance);
fprintf('Theoretical acceptance fraction: %.4f\n', theoreticalAcceptance);
fprintf('Mean time per accepted sample: %.8e seconds\n', meanTime);
fprintf('Histogram bin width: %.2f\n', binWidth);

histogram(samples, ...
    'BinWidth', binWidth, ...
    'Normalization', 'pdf');

hold on

x = linspace(0, max(samples), 1000);
f = x .* exp(-x);

plot(x, f, 'LineWidth', 2)

xlabel('x')
ylabel('Probability Density')
title(sprintf('Gamma(2,1) Rejection Sampling, N = %d', N))
legend('Empirical PDF', 'Theoretical f(x) = xe^{-x}')

hold off