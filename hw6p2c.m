clear; clc; close all;

D = readtable('big12.csv');

D.home = string(D.home);
D.away = string(D.away);

teams = unique([D.home; D.away]);
m = numel(teams);
w = 1;
targetWeek = 6;

[~, home] = ismember(D.home, teams);
[~, away] = ismember(D.away, teams);

h = @(x) 0.5 + 0.5.*sign(x-0.5).*sqrt(abs(2*x-1));

H25 = zeros(m);
C25 = zeros(m);
wins = zeros(m,1);
games = zeros(m,1);

for t = find(D.season == 2025 & D.completed == 1)'
    i = home(t);
    j = away(t);
    S = D.home_pts(t);
    T = D.away_pts(t);

    v = h((S+1)/(S+T+2));

    H25(i,j) = H25(i,j) + v;
    H25(j,i) = H25(j,i) + 1-v;
    C25(i,j) = C25(i,j) + 1;
    C25(j,i) = C25(j,i) + 1;

    games([i j]) = games([i j]) + 1;
    wins(i) = wins(i) + (S>T) + 0.5*(S==T);
    wins(j) = wins(j) + (T>S) + 0.5*(S==T);
end

flat = 0.5*ones(m);

lastSeason = flat;
met = C25 > 0;
lastSeason(met) = H25(met)./C25(met);

s = (wins+1)./(games+2);
custom = s./(s+s');

priors = cat(3, flat, lastSeason, custom);
names = ["Flat"; "2025 head-to-head"; "2025 win strength"];

H26 = zeros(m);
n = zeros(m,1);

for t = find(D.season == 2026 & ...
             D.completed == 1 & D.week < targetWeek)'
    i = home(t);
    j = away(t);
    S = D.home_pts(t);
    T = D.away_pts(t);

    v = h((S+1)/(S+T+2));

    H26(i,j) = H26(i,j) + v;
    H26(j,i) = H26(j,i) + 1-v;
    n([i j]) = n([i j]) + 1;
end

R = zeros(m,3);
lambdaStar = zeros(3,1);
ratio = zeros(3,1);

for a = 1:3
    Pi = priors(:,:,a);
    Pi(1:m+1:end) = 0;

    A = (H26+w*Pi)./(n+w*(m-1));

    [V,L] = eig(A);
    z = diag(L);
    [~, order] = sort(abs(z), 'descend');

    lambdaStar(a) = real(z(order(1)));
    ratio(a) = abs(z(order(2)))/lambdaStar(a);

    r = real(V(:,order(1)));
    if sum(r) < 0
        r = -r;
    end
    R(:,a) = r/sum(r);
end

disp(table(names, lambdaStar, ratio, ...
    'VariableNames', ...
    {'Prior','LambdaStar','AbsLambda2OverLambdaStar'}));

ranking = table(teams, R(:,1), R(:,2), R(:,3), ...
    'VariableNames', {'Team','Flat','LastSeason','Custom'});
disp(sortrows(ranking, 'Flat', 'descend'));

[~, order] = sort(R(:,1), 'descend');

figure;
bar(R(order,:), 'grouped');
xticks(1:m);
xticklabels(teams(order));
xtickangle(60);
ylabel('Normalized strength r_i');
title('Big 12 rankings before Oct 10-11, 2026 (w = 1)');
legend(names, 'Location', 'best');
grid on;

rows = find(D.season == 2026 & D.week == targetWeek);

predictions = table(D.home(rows), D.away(rows), ...
    'VariableNames', {'Home','Away'});

picks = strings(numel(rows),3);

for a = 1:3
    difference = R(home(rows),a) - R(away(rows),a);

    picks(:,a) = D.home(rows);
    picks(difference < 0,a) = D.away(rows(difference < 0));
    picks(abs(difference) < 1e-12,a) = "Tie";
end

predictions.Flat = picks(:,1);
predictions.LastSeason = picks(:,2);
predictions.Custom = picks(:,3);

actual = repmat("Pending", numel(rows), 1);

for g = 1:numel(rows)
    t = rows(g);

    if D.completed(t) == 1
        if D.home_pts(t) > D.away_pts(t)
            actual(g) = D.home(t);
        elseif D.home_pts(t) < D.away_pts(t)
            actual(g) = D.away(t);
        else
            actual(g) = "Tie";
        end
    end
end

predictions.Actual = actual;
disp(predictions);

done = actual ~= "Pending";

if any(done)
    for a = 1:3
        correct = sum(picks(done,a) == actual(done));
        fprintf('%s: %d correct out of %d completed games\n', ...
            char(names(a)), correct, sum(done));
    end
else
    disp('Week 6 results are still unplayed in this CSV.');
end