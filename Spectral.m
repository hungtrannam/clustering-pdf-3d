clear; close all;
numCluster = 3;
h = .01;
x = -5: h : 15;

mu = [0.3 4.0 9.1 1.0 5.5 8.0 4.8];
numSample = length(mu);
f = [];
for i = 1:length(mu)
    f = [f normpdf(x, mu(i), 1)'];
end

sigma = 1; 
for j = 1:numSample
    for i = 1:numSample
        Wf(i, j) = Integration(h, abs(f(:, i) - f(:, j)), 1) + 10^(-10);
    end
end

similarityMatrix = exp(-Wf.^2 / (2 * sigma^2));

D = diag(sum(similarityMatrix, 2));
L = D - similarityMatrix;
L_norm = inv(sqrt(D)) * L * inv(sqrt(D));

[eigenVectors, ~] = eigs(L_norm, numCluster, 'smallestabs');

[idx, ~] = kmeans(eigenVectors, numCluster);

color = {'#00AFEF', '#1086CC', '#1F5CA9'};
for i = 1:numCluster
    plot(x, f(:, idx == i), 'Color', color{i}); hold on;
    plot(x, mean(f(:,idx == i), 2), 'Color', color{i}, 'LineWidth', 2.5);
    hold on;
end
box off;
set(gca, 'YColor', 'none');
