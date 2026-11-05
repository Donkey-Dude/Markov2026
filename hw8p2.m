n = 2e5;

s1 = binornd(poissrnd(15,n,1),0.4);
b1 = binornd(poissrnd(10,n,1),0.4);

s2 = poissrnd(6,n,1);
b2 = poissrnd(4,n,1);
total = s2+b2;

disp([mean(s1) var(s1); mean(b1) var(b1); mean(s2) var(s2); mean(b2) var(b2)])
disp(corr(s2,b2))

figure
histogram(total,'Normalization','probability')
hold on
k = 0:25;
plot(k,poisspdf(k,10),'ro-')

disp(mean(s2(total==10)==6))

ts = -log(rand(n,1))/1.2;
tb = -log(rand(n,1))/0.8;

first = min(ts,tb);
starFirst = ts<tb;
M = max(ts,tb);

disp(mean(starFirst))
disp(mean(starFirst & first<=1))
disp(mean(M))
disp(corr(first,double(starFirst)))