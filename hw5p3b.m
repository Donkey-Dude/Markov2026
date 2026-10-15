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

q = ones(1,8)/8;
iterations = 0;

while true
    qNext = q*G;
    change = norm(qNext-q,1);
    q = qNext;
    iterations = iterations+1;

    if change < 1e-10
        break
    end
end

A = G' - eye(8);
A(end,:) = ones(1,8);
rhs = [zeros(7,1); 1];
pi = (A \ rhs)';

fprintf('Iterations: %d\n', iterations);
fprintf('Final change: %.12e\n', change);
fprintf('Difference from linear solve: %.12e\n', norm(q-pi,1));
format long g
disp('Stationary distribution:')
disp(pi)