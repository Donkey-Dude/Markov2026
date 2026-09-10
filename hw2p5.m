N = 1e5;

U = rand(N,1);
Z = 1 - sqrt(1-U);

sampleMean = mean(Z);

fprintf('Sample mean: %.6f\n', sampleMean);