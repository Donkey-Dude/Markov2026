R = 2e4;
T = 1e4;
d0 = 10;

lamb = zeros(R,1);
lions = d0 * ones(R,1);
alive = true(R,1);

S1 = zeros(1,T);

for t = 1:T
    lambMove = (rand(R,1) < 0.5) * 2 - 1;
    lamb(alive) = lamb(alive) + lambMove(alive);

    lionMove = (rand(R,1) < 0.5) * 2 - 1;
    lions(alive) = lions(alive) + lionMove(alive);

    aliveIdx = find(alive);
    caught = lions(alive) == lamb(alive);
    alive(aliveIdx(caught)) = false;

    S1(t) = sum(alive)/R;
end

time = 1:T;

valid1 = (time >= 1e2) & (time <= 1e4) & (S1 > 0);

x1 = log(time(valid1));
y1 = log(S1(valid1));

p1 = polyfit(x1, y1, 1);
beta1 = -p1(1);

yfit1 = p1(1)*x1 + p1(2);
Sfit1 = exp(yfit1);

positive1 = S1 > 0;

S_theory = erf(d0./(2*sqrt(time)));

lamb = zeros(R,1);
lions = d0 * ones(R,2);
alive = true(R,1);

S2 = zeros(1,T);

for t = 1:T
    lambMove = (rand(R,1) < 0.5) * 2 - 1;
    lamb(alive) = lamb(alive) + lambMove(alive);

    lionMove = (rand(R,2) < 0.5) * 2 - 1;
    lions(alive,:) = lions(alive,:) + lionMove(alive,:);

    aliveIdx = find(alive);
    caught = any(lions(alive,:) == lamb(alive), 2);
    alive(aliveIdx(caught)) = false;

    S2(t) = sum(alive)/R;
end

valid2 = (time >= 1e2) & (time <= 1e4) & (S2 > 0);

x2 = log(time(valid2));
y2 = log(S2(valid2));

p2 = polyfit(x2, y2, 1);
beta2 = -p2(1);

positive2 = S2 > 0;

S1_squared = S1.^2;

loglog(time(positive1), S1(positive1))
hold on

loglog(time, S_theory)
loglog(time(valid1), Sfit1)
loglog(time(positive2), S2(positive2))
loglog(time, S1_squared)

xlabel('Time t')
ylabel('Survival probability')

title(sprintf('beta_1 = %.3f, beta_2 = %.3f, fit: 10^2 <= t <= 10^4', ...
    beta1, beta2))

legend('S_1 simulation', ...
       'S_1 continuum prediction', ...
       'S_1 power-law fit', ...
       'S_2 simulation', ...
       'S_1^2')

hold off

tableTimes = [100 1000 10000];

S2_values = S2(tableTimes);
S1_squared_values = S1_squared(tableTimes);

results = table(tableTimes', S2_values', S1_squared_values', ...
    'VariableNames', {'Time', 'S2', 'S1_squared'});

disp(results)