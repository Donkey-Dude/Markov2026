
rng(1)

n = [10 100 1000];
meanScatter = zeros(1,3);
medianScatter = zeros(1,3);

for i = 1:3
    x = 1.5 + 0.8*tan(pi*(rand(n(i),2000)-0.5));

    m = mean(x,1);
    d = median(x,1);

    meanScatter(i) = (prctile(m,75)-prctile(m,25))/1.349;
    medianScatter(i) = (prctile(d,75)-prctile(d,25))/1.349;
end

disp(table(n',meanScatter',medianScatter', ...
    'VariableNames',{'n','MeanScatter','MedianScatter'}))