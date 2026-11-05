T = 720;
lambda = @(t) 12 + 8*cos(pi*(t-16)/12);

n = poissrnd(20*T);
proposals = sort(T*rand(n,1));
times = proposals(rand(n,1) < lambda(proposals)/20);

disp(length(times)/n)

hours = mod(times,24);
counts = histcounts(hours,0:24)/30;

figure
bar(0.5:23.5,counts)
hold on
t = linspace(0,24,500);
plot(t,lambda(t),'r','LineWidth',2)
xlabel('Hour')
ylabel('Arrivals per hour per day')
legend('Simulated','True rate')

days = floor(times/24)+1;
early = accumarray(days(hours>=3 & hours<4),1,[30 1]);

disp([mean(early) var(early)])
disp([4.091 4.091])

logL = @(p) sum(log(p(1)+p(2)*cos(pi*(times-p(3))/12))) - T*p(1);

objective = @(p) -logL(p) + 1e10*(p(1)<=abs(p(2)));

p = fminsearch(objective,[10 5 14]);

Lfit = logL(p);
r = length(times)/T;
Lconst = length(times)*log(r)-r*T;

disp(p)
disp([Lfit Lconst])