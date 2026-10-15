P = zeros(5);
for k = 0:4
    if k < 4
        P(k+1, k+2) = (4-k)/4;
    end
    if k > 0
        P(k+1, k) = k/4;
    end
end

Q = zeros(61,5);
Q(1,:) = [1 0 0 0 0];

for n = 1:60
    Q(n+1,:) = Q(n,:) * P;
end

format long
disp('q_50:')
disp(Q(51,:))
disp('q_51:')
disp(Q(52,:))

n = (0:60)';
q2 = Q(:,3);
q4 = Q(:,5);
avg2 = cumsum(q2) ./ (n+1);

figure;
plot(n, q2, '-o', n, q4, '-s', n, avg2, '-', ...
    'LineWidth', 1.2, 'MarkerSize', 3);
hold on;
yline(3/8, '--', '\pi(2) = 3/8');
xlabel('n');
ylabel('Probability / running average');
legend('q_n(2)', 'q_n(4)', ...
    'Running average of q_k(2)', '\pi(2)', ...
    'Location', 'best');
grid on;