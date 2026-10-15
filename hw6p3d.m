clear; clc;
rng(12345);

N = 100;
K = 10;
gamma = 1;
dt = 0.005;
Rvalues = [0.5; 0.8; 1.5; 3];

n = (1:K-1)';
expectedTime = zeros(size(Rvalues));

for k = 1:numel(Rvalues)
    beta = Rvalues(k)*gamma/N;

    p = beta*n.*(N-n)*dt;
    q = gamma*n*dt;

    A = diag(p+q) ...
        - diag(p(1:end-1), 1) ...
        - diag(q(2:end), -1);

    e = A \ ones(K-1,1);
    expectedTime(k) = e(1)*dt;
end

disp(table(Rvalues, expectedTime, ...
    'VariableNames', {'R0','ExpectedPhysicalTime'}));

R0 = 1.5;
beta = R0*gamma/N;
numTrials = 100000;

X = ones(numTrials,1);
steps = zeros(numTrials,1);

while true
    active = find(X > 0 & X < K);
    if isempty(active)
        break
    end

    state = X(active);
    p = beta*state.*(N-state)*dt;
    q = gamma*state*dt;
    u = rand(numel(active),1);

    X(active) = state + double(u < p) - double(u >= 1-q);
    steps(active) = steps(active) + 1;
end

times = steps*dt;
estimate = mean(times);
SE = std(times)/sqrt(numTrials);

fprintf('Exact at R0 = 1.5: %.6f\n', expectedTime(3));
fprintf('Simulation: %.6f, standard error %.6f\n', estimate, SE);
fprintf('Approximate 95%% CI: [%.6f, %.6f]\n', ...
    estimate-1.96*SE, estimate+1.96*SE);