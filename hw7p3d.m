
rng(1)

lambda = 4;
R = 200000;

N1 = zeros(R,1);
N2 = poissrnd(lambda,R,1);
T3 = zeros(R,1);
first1 = zeros(R,1);
half1 = zeros(R,1);
half2 = zeros(R,1);
halfU = zeros(R,1);
firstU = nan(R,1);

figure
hold on

for i = 1:R
    times = cumsum(-log(rand(1,30))/lambda);
    T3(i) = times(3);
    arrivals = times(times<=1);
    N1(i) = length(arrivals);
    half1(i) = sum(arrivals<=0.5);
    half2(i) = sum(arrivals>0.5);

    if N1(i)>0
        first1(i) = arrivals(1);
    end

    u = sort(rand(1,N2(i)));
    halfU(i) = sum(u<=0.5);

    if N2(i)>0
        firstU(i) = u(1);
    end

    if i<=10
        subplot(2,1,1)
        hold on
        plot(arrivals,i*ones(size(arrivals)),'k|')
        subplot(2,1,2)
        hold on
        plot(u,i*ones(size(u)),'k|')
    end
end

subplot(2,1,1)
title('Exponential gaps')
xlim([0 1])
ylim([0 11])
subplot(2,1,2)
title('Poisson + uniform')
xlim([0 1])
ylim([0 11])

k = 0:max([N1;N2]);
pmf = exp(-4)*4.^k./factorial(k);

figure
histogram(N1,'BinEdges',-0.5:max(k)+0.5,'Normalization','probability')
hold on
plot(k,pmf,'ro-')
title('Exponential gaps')
legend('Simulation','Poisson(4)')

figure
histogram(N2,'BinEdges',-0.5:max(k)+0.5,'Normalization','probability')
hold on
plot(k,pmf,'ro-')
title('Poisson + uniform')
legend('Simulation','Poisson(4)')

disp(table([mean(N1);mean(N2)],[var(N1);var(N2)], ...
    'VariableNames',{'Mean','Variance'}, ...
    'RowNames',{'Gaps','Uniform'}))

t = linspace(0,max(T3),500);
f = lambda^3*t.^2.*exp(-lambda*t)/2;

figure
histogram(T3,80,'Normalization','pdf')
hold on
plot(t,f,'LineWidth',2)
title('Third arrival time')
legend('Simulation','Gamma density')

idx1 = N1==3;
idx2 = N2==3;

disp(['Gaps P(N(0.5)=0 | N(1)=3): ' num2str(mean(half1(idx1)==0))])
disp(['Uniform P(N(0.5)=0 | N(1)=3): ' num2str(mean(halfU(idx2)==0))])
disp(['Gaps E[T1 | N(1)=3]: ' num2str(mean(first1(idx1)))])
disp(['Uniform E[T1 | N(1)=3]: ' num2str(mean(firstU(idx2)))])
disp(['Half-hour correlation: ' num2str(corr(half1,half2))])

m = round(logspace(log10(10),log10(2000),40));
m = unique(m);
dt = 1./m;
TV = zeros(size(m));

for i = 1:length(m)
    n = m(i);
    j = 0:n;
    p = exp(gammaln(n+1)-gammaln(j+1)-gammaln(n-j+1) ...
        +j*log(lambda/n)+(n-j)*log(1-lambda/n));
    q = exp(-lambda+j*log(lambda)-gammaln(j+1));
    TV(i) = 0.5*(sum(abs(p-q))+1-sum(q));
end

c = polyfit(log(dt),log(TV),1);

figure
loglog(dt,TV,'o')
hold on
loglog(dt,exp(polyval(c,log(dt))),'LineWidth',2)
xlabel('Delta t')
ylabel('Total variation distance')
legend('TV distance','Power-law fit')

disp(['Log-log slope: ' num2str(c(1))])
