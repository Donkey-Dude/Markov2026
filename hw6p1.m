ms = (3:12)';
M = NaN(size(ms));
lambda2 = zeros(size(ms));
bound = (ms - 1).^2 + 1;

for k = 1:numel(ms)
    m = ms(k);

    P = zeros(m);
    for i = 1:m-1
        P(i,i+1) = 1;
    end
    P(m,1) = 1/2;
    P(m,3) = 1/2;

    B = eye(m);
    for n = 1:bound(k)
        B = double(B * P > 0);
        if all(B(:))
            M(k) = n;
            break
        end
    end

    moduli = sort(abs(eig(P)), 'descend');
    lambda2(k) = moduli(2);
end

results = table(ms, M, bound, lambda2, ...
    'VariableNames', {'m','M','WielandtBound','AbsLambda2'});
disp(results)