
rng(1)

lambda = [3 2 1.5 1.2 1 0.75];
n = length(lambda);
N = 200000;

E = zeros(2^n,1);
P = zeros(2^n,n);

for mask = 1:2^n-1
    ids = find(bitget(mask,1:n));
    r = sum(lambda(ids));
    E(mask+1) = 1/r;

    if length(ids)==1
        P(mask+1,ids) = 1;
    else
        for i = ids
            next = bitset(mask,i,0);
            E(mask+1) = E(mask+1) + lambda(i)/r*E(next+1);
            P(mask+1,:) = P(mask+1,:) + lambda(i)/r*P(next+1,:);
        end
    end
end

exactE = E(end);
lastProb = P(end,:);

altE = 0;
for mask = 1:2^n-1
    ids = find(bitget(mask,1:n));
    altE = altE + (-1)^(length(ids)+1)/sum(lambda(ids));
end

X = -log(rand(N,n))./lambda;
M = max(X,[],2);
[~,last] = max(X,[],2);
[firstTime,first] = min(X,[],2);

simProb = histcounts(last,0.5:1:n+0.5)/N;

disp(table((1:n)',lastProb',simProb', ...
    'VariableNames',{'Student','ExactLast','SimulatedLast'}))

disp(table(exactE,altE,mean(M), ...
    'VariableNames',{'GroupMethod','Alternating','Simulation'}))

disp(table((1:n)',accumarray(first,firstTime,[n 1],@mean), ...
    'VariableNames',{'FirstStudent','MeanFirstTime'}))

t = linspace(0,max(M),500);
F = prod(1-exp(-t(:)*lambda),2);
f = zeros(size(t(:)));

for i = 1:n
    others = setdiff(1:n,i);
    f = f + lambda(i)*exp(-lambda(i)*t(:)) ...
        .*prod(1-exp(-t(:)*lambda(others)),2);
end

figure
histogram(M,80,'Normalization','pdf')
hold on
plot(t,f,'LineWidth',2)
xlabel('Last finishing time')
ylabel('Density')
legend('Simulation','Exact')

K = 12;
EK = zeros(1,K);

for k = 1:K
    rates = 1:k;
    for mask = 1:2^k-1
        ids = find(bitget(mask,1:k));
        EK(k) = EK(k) + (-1)^(length(ids)+1)/sum(rates(ids));
    end
end

figure
plot(1:K,EK,'o-')
hold on
yline(2*log(2),'--')
xlabel('Number of students')
ylabel('E[M]')
legend('Exact E[M]','2 ln 2')
