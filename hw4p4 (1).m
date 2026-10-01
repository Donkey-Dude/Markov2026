clear; clc; close all;
rng(1);

P = [0   1/2 1/2 0   0;
    1/4 0   0   1/2 1/4;
    3/4 0   0   0   1/4;
    0   0   0   1   0;
    0   0   0   0   1];

Q = P(1:3,1:3);
R = P(1:3,4:5);
CDF = cumsum(P,2);

exact = [1/2 4 4   4;
    5/8 2 9/5 7/3;
    3/8 4 5   17/5];

T = zeros(10^4,3);
fate = zeros(10^4,3);
estimate = zeros(3,4);

for x = 1:3
    for j = 1:10^4
        s = x;
        while s <= 3
            s = find(rand < CDF(s,:),1);
            T(j,x) = T(j,x) + 1;
        end
        fate(j,x) = s;
    end

    folded = fate(:,x) == 4;
    estimate(x,:) = [mean(folded), mean(T(:,x)), ...
        mean(T(folded,x)), mean(T(~folded,x))];
end

disp(array2table([exact estimate], ...
    'RowNames',{'U','I','M'}, ...
    'VariableNames',{'h','g','tauF','tauA', ...
    'h_sim','g_sim','tauF_sim','tauA_sim'}));

n = 1:max(T(:,2));
pmf = zeros(length(n),2);

for k = n
    joint = Q^(k-1)*R;
    pmf(k,:) = joint(2,:) ./ [5/8 3/8];
end

figure;
labels = {'Folded','Aggregated'};

for f = 1:2
    subplot(1,2,f);
    times = T(fate(:,2) == f+3,2);

    histogram(times,0.5:1:(n(end)+0.5), ...
        'Normalization','probability');
    hold on;
    stem(n,pmf(:,f),'k.','LineWidth',1);

    xlabel('Absorption time T');
    ylabel('Conditional probability');
    title(labels{f});
    legend('Simulation','Exact PMF');
    grid on;
end