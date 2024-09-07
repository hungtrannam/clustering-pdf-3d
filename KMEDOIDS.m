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
initial_idx = randperm(numSample, numCluster);
centroids = f(:, initial_idx);

distances = zeros(numCluster, numSample);
while iter < max_iter
    iter = iter + 1;
    
    % Compute distances
    for j = 1:numSample
        for i = 1:numCluster
            distances(i, j) = Integration(h, abs(centroids(:, i) - f(:, j)), 1) + 10^(-10);
        end
    end
    [~, IDX] = min(distances);
    
    % K-medoids update step
    new_centroids_idx = zeros(1, numCluster);
    for i = 1:numCluster
        cluster_points = find(IDX == i);
        if ~isempty(cluster_points)
            for u = 1:numCluster-1
                for v = 1:numCluster-1
                    pairwise_dist(u, v) = Integration(h, abs(f(:,cluster_points(u)) - f(:,cluster_points(v))), 1);
                end
            end
            total_dist = sum(pairwise_dist);
            [~, min_idx] = min(total_dist);
            new_centroids_idx(i) = cluster_points(min_idx);
        end
    end
    
    centroidsnew = f(:, new_centroids_idx);
    
    ObjFun = sum(min(distances));
    fprintf('Iteration count = %d, obj. kmedoids = %f\n', iter, ObjFun);
    
    if iter > 1 && norm(centroidsnew - centroids,1) < epsilon
        break;
    end
    centroids = centroidsnew;

    for i = 1:numCluster
        plot(x, f(:,IDX == i),'Color', color{i}, 'LineStyle',':', 'LineWidth', 2.5, 'DisplayName', 'Density Function Vectors'); hold on;
        plot(x, centroids(:,i), 'LineWidth', 3, 'Color', color{i}, 'DisplayName', 'Centroids');
        box off;
        set(gca, 'YColor', 'none');
    end
end
