R = 2e4;
T = 1e4;
d0 = 10;
lamb = zeros(R,1);
lions = 10 * ones(R,1);
alive = true(R,1);
S = zeros(1,T);

for t = 1:T
    lambMove = (rand(R,1) < 0.5) * 2 - 1;
    lamb(alive) = lamb(alive) + lambMove(alive);

    lionMove = (rand(R,1) < 0.5) * 2 - 1;
    lions(alive) = lions(alive) + lionMove(alive);

    aliveIdx = find(alive);
    caught = lions(alive) == lamb(alive);
    alive(aliveIdx(caught)) = false;
    S(t) = sum(alive)/R;
end

time = 1:T;

valid = (time >= 1e2) & (time <= 1e4) & (S > 0);

x = log(time(valid));
y = log(S(valid));

p = polyfit(x, y, 1);
beta = -p(1);

positive = S > 0;

S_theory = erf(d0./(2*sqrt(time)));

yfit = p(1)*x + p(2);
Sfit = exp(yfit);

loglog(time(positive), S(positive))
hold on

loglog(time, S_theory)
loglog(time(valid), Sfit)

xlabel('Time t')
ylabel('Survival probability S_1(t)')
title(sprintf('N = 1, beta = %.3f, fit: 10^2 <= t <= 10^4', beta))
legend('Simulation', 'Continuum prediction', 'Power-law fit')

hold off