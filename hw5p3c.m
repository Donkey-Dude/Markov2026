clear; clc;
rng(1);                       

P = [
    0   1/2 1/2 0   0   0   0   0;
    0   0   1/2 1/2 0   0   0   0;
    1/2 0   0   0   1/2 0   0   0;
    1/3 0   1/3 0   1/3 0   0   0;
    0   1/3 0   0   0   1/3 1/3 0;
    0   0   0   0   0   1   0   0;
    0   0   0   0   0   0   0   1;
    0   0   0   0   0   0   1   0
    ];

d = 0.85;
G = d*P + (1-d)/8 * ones(8);

A = G' - eye(8);
A(end,:) = ones(1,8);
pi = (A \ [zeros(7,1); 1])';

CDF = cumsum(G,2);
CDF(:,end) = 1;

R = 100;
Tmax = 1e5;
state = ones(R,1);             
counts = zeros(R,8);
surfer = (1:R)';

times = unique(round(logspace(0,5,101)));
rmsError = zeros(size(times));
checkpoint = 1;

for t = 1:Tmax
    u = rand(R,1);
    state = 1 + sum(u > CDF(state,:), 2);

    idx = sub2ind([R,8], surfer, state);
    counts(idx) = counts(idx) + 1;

    if t == times(checkpoint)
        estimates = counts/t;
        maxError = max(abs(estimates-pi), [], 2);
        rmsError(checkpoint) = sqrt(mean(maxError.^2));
        checkpoint = checkpoint + 1;
    end
end

piHat = counts/Tmax;

fitMask = times >= 1e2 & times <= 1e5;
coeff = polyfit(log(times(fitMask)), log(rmsError(fitMask)), 1);
slope = coeff(1);
fittedError = exp(polyval(coeff, log(times(fitMask))));

iidBenchmark = sqrt(pi(6)*(1-pi(6))/Tmax);

fprintf('Fitted slope: %.4f\n', slope);
fprintf('RMS maximum error at T = 10^5: %.6g\n', rmsError(end));
fprintf('Independent-sample benchmark: %.6g\n', iidBenchmark);
fprintf('Ratio to benchmark: %.3f\n', rmsError(end)/iidBenchmark);

figure;
tiledlayout(1,2);

nexttile;
bar(1:8, [pi(:), piHat(1,:)']);
xlabel('Page');
ylabel('Probability');
title('One surfer at T = 10^5');
legend('\pi', 'Time-average estimate', 'Location', 'best');
grid on;

nexttile;
loglog(times, rmsError, 'o-', 'MarkerSize', 3);
hold on;
loglog(times(fitMask), fittedError, '--', 'LineWidth', 1.5);
xlabel('T');
ylabel('RMS maximum error');
title(sprintf('Fitted slope = %.3f', slope));
legend('Simulation', 'Fit over 10^2 to 10^5', ...
    'Location', 'best');
grid on;