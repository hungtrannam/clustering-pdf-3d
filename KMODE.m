clear; close all;
numCluster = 3;
max_iter = 10; 
iter = 0;
epsilon = 0.01; 
color = {'#00AFEF', '#1086CC', '#1F5CA9'};

h = .01;
x = -5: h : 15;

mu = [0.3 4.0 9.1 1.0 5.5 8.0 4.8];
numSample = length(mu);

f = [];
for i = 1:length(mu)
    f = [f normpdf(x, mu(i), 1)'];
end

rng(2);
centroids = f(:, randperm(numSample, numCluster));

while iter < max_iter
    iter = iter + 1;
    
    for j = 1:numSample
        for i = 1:numCluster
            distances(i, j) = Integration(h, abs(centroids(:, i) - f(:, j)), 1) + 10^(-10);
        end
    end
    [~, IDX] = min(distances);

    for i = 1:numCluster
        if sum(IDX == i) > 0
            centroids(:, i) = mode(f(:, IDX == i),2);
        end
    end

    ObjFun = sum(min(distances));
    fprintf('Iteration count = %d, obj. kmeans = %f\n', iter, ObjFun);

    if iter > 1 && norm(centroidsnew - centroids) < epsilon
        break;
    end
    centroidsnew = centroids;

    % Hiển thị kết quả
    for i = 1:numCluster
        plot(x, f(:,IDX == i),'Color', color{i}, 'LineStyle',':', 'LineWidth', 2.5, 'DisplayName', 'Density Function Vectors'); hold on;
        plot(x, centroids(:,i), 'LineWidth', 3, 'Color', color{i}, 'DisplayName', 'Centroids');
        box off;
        set(gca, 'YColor', 'none');
    end
end
