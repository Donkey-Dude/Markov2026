clear; clc;

a = 0.3;
b = 0.1;

P = zeros(5);
for k = 0:4
    if k < 4
        P(k+1,k+2) = a*(4-k)/4;
    end
    if k > 0
        P(k+1,k) = b*k/4;
    end
    P(k+1,k+1) = 1 - sum(P(k+1,:));
end

A = P' - eye(5);
A(end,:) = ones(1,5);
rhs = [0; 0; 0; 0; 1];
pi = (A \ rhs)';

theta = a/(a+b);
piBin = zeros(1,5);
for k = 0:4
    piBin(k+1) = nchoosek(4,k)*theta^k*(1-theta)^(4-k);
end

disp('Stationary distribution:')
disp(pi)
fprintf('Difference from binomial formula: %.3e\n', ...
    max(abs(pi-piBin)));

Q = [1 0 0 0 0];
err = max(abs(Q-pi));
n = 0;

while err(end) >= 1e-6
    Q(end+1,:) = Q(end,:)*P;
    n = n+1;
    err(end+1,1) = max(abs(Q(end,:)-pi));
end

fprintf('Smallest n: %d\n', n);
fprintf('Error at n-1: %.12e\n', err(end-1));
fprintf('Error at n:   %.12e\n', err(end));

moduli = sort(abs(eig(P)), 'descend');
fprintf('Second largest eigenvalue modulus: %.12f\n', moduli(2));

figure;
h = plot(0:n, Q, 'LineWidth', 1.4);
hold on;
for k = 1:5
    yline(pi(k), '--', 'Color', h(k).Color, ...
        'HandleVisibility', 'off');
end
xlabel('n');
ylabel('Probability');
title('a = 0.3, b = 0.1');
legend('q_n(0)', 'q_n(1)', 'q_n(2)', 'q_n(3)', 'q_n(4)', ...
    'Location', 'best');
grid on;