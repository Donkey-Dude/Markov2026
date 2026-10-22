
rng(1)

x = [-1.77 2.47 1.34 2.18 13.05 1.60 1.50 -1.97 0.79 1.50 2.00 2.63 1.19 -2.31 0.87];
n = length(x);
T = 50000;
svals = [0.02 0.5 5];

logg = @(a,b) n*b - sum(log(exp(2*b) + (x-a).^2));

chains = cell(3,1);
rates = zeros(3,1);
tau = zeros(3,1);

figure
for j = 1:3
    s = svals(j);
    z = zeros(T,2);
    accepted = 0;

    for i = 2:T
        y = z(i-1,:) + s*randn(1,2);
        if log(rand) < logg(y(1),y(2)) - logg(z(i-1,1),z(i-1,2))
            z(i,:) = y;
            accepted = accepted + 1;
        else
            z(i,:) = z(i-1,:);
        end
    end

    chains{j} = z;
    rates(j) = accepted/(T-1);

    a = z(5001:end,1);
    a = a - mean(a);
    c = real(ifft(abs(fft(a,2^nextpow2(2*length(a)))).^2));
    c = c(1:length(a))./(length(a):-1:1)';
    c = c/c(1);
    k = find(c(2:end)<=0,1);
    if isempty(k), k = length(c)-1; end
    tau(j) = 1 + 2*sum(c(2:k));

    subplot(3,1,j)
    plot(z(1:400,1))
    title(['s = ' num2str(s)])
    ylabel('alpha')
end
xlabel('Step')

disp(table(svals',rates,tau,'VariableNames',{'s','Acceptance','IACT'}))

[~,best] = min(tau);
z = chains{best}(5001:end,:);
alpha = z(:,1);
beta = exp(z(:,2));

A = linspace(-5,8,300);
B = linspace(-4,3,300);
[AA,BB] = meshgrid(A,B);
L = zeros(size(AA));

for i = 1:numel(AA)
    L(i) = logg(AA(i),BB(i));
end

P = exp(L-max(L(:)));
P = P/sum(P(:));

figure
contour(AA,BB,P,20)
hold on
plot(z(1:10:end,1),z(1:10:end,2),'.','MarkerSize',2)
xlabel('alpha')
ylabel('log beta')
title(['s = ' num2str(svals(best))])

gridAlpha = AA(:);
gridBeta = exp(BB(:));
weights = P(:);

[sa,ia] = sort(gridAlpha);
ca = cumsum(weights(ia));
[ua,idxa] = unique(ca);
qa = interp1(ua,sa(idxa),[0.05 0.95]);

[sb,ib] = sort(gridBeta);
cb = cumsum(weights(ib));
[ub,idxb] = unique(cb);
qb = interp1(ub,sb(idxb),[0.05 0.95]);

disp(['Best s: ' num2str(svals(best))])
disp(['Chain alpha: ' num2str(mean(alpha)) '  ' mat2str(prctile(alpha,[5 95]))])
disp(['Grid alpha:  ' num2str(sum(gridAlpha.*weights)) '  ' mat2str(qa)])
disp(['Chain beta:  ' num2str(mean(beta)) '  ' mat2str(prctile(beta,[5 95]))])
disp(['Grid beta:   ' num2str(sum(gridBeta.*weights)) '  ' mat2str(qb)])
disp(['Sample mean:   ' num2str(mean(x))])
disp(['Sample median: ' num2str(median(x))])
