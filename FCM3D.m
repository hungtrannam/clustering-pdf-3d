clear; close all;
m = 2;
numCluster = 2;
max_iter = 100;
iter = 0;
epsilon = 0.001;
color2 = {'#00AFEF', '#1F5CA9'};
color3 = {'#00AFEF', '#1086CC', '#1F5CA9'};
color = {'#00AFEF', '#05A1E3', '#0A93D8', '#1086CC', '#1578C0', '#1A6AB5', '#1F5CA9'};

h = 0.05;
[x, y] = meshgrid(-4:h:11, -4:h:11);

mu = [0 0; 2 0; 1 sqrt(3); 5 7; 7 7; 5 5; 7 5];

rng(1);

fi = [];
figure(1);
hold on;
for i = 1:size(mu, 1)
    fi(:,:,i) = mvnpdf([x(:) y(:)], mu(i,:), eye(2));
    f(:,:,i) = reshape(fi(:,:,i), size(x,1), size(x,2));
    contour(x, y, f(:,:,i)); hold on;
end
hold off;

numSample = size(mu, 1);

U0 = rand(numCluster, numSample);
U0 = U0 ./ sum(U0);

% Initial fv calculation
    % Update fvnew
    for i = 1:numCluster
        weights = U0(i, :).^m;
        normalization_factor = sum(weights);
        weighted_sum = sum(f .* reshape(weights, 1, 1, numSample), 3);
        fv(:,:,i) = weighted_sum / normalization_factor;
    end
fv_copy = fv;

while iter < max_iter
    iter = iter + 1;

    % Compute Wf
    Wf = zeros(numCluster, numSample);
    for j = 1:numSample
        for i = 1:numCluster
            diff = abs(fv(:,:,i) - f(:,:,j));
            Wf(i, j) = (Integration(h, diff, 2)) + 10^(-10);
        end
    end

    % Update Unew
    Unew = zeros(numCluster, numSample);
    for j = 1:numSample
        for i = 1:numCluster
            numerator = 1 / (Wf(i, j)) ^ (2/(m - 1));
            denominator = 0;
            for k = 1:numCluster
                denominator = denominator + 1 / (Wf(k, j)) ^ (2/(m - 1));
            end
            Unew(i, j) = numerator / denominator;
        end
    end

    % Update fvnew
    for i = 1:numCluster
        weights = Unew(i, :).^m;
        normalization_factor = sum(weights);
        weighted_sum = sum(f .* reshape(weights, 1, 1, numSample), 3);
        fvnew(:,:,i) = weighted_sum / normalization_factor;
    end

    % Compute objective function value
    ObjFun = sum(sum((Unew) .* Wf));
    fprintf('Iteration count = %d, obj. fcm = %f\n', iter, ObjFun);


    % Check convergence
    if norm(U0 - Unew, 1) < epsilon
        break;
    end

    % Update for the next iteration
    fv = fvnew;
    U0 = Unew;


end

% Plotting results
figure(3);
for i = 1:numCluster
    meshc(x, y, fvnew(:,:,i)); hold on;
end

figure(6);
for j = 1:numSample
    meshc(x, y, f(:,:,j)); hold on;
end

hold off;

figure(4);
heatmap(Unew);


figure(5);
for i = 1:numCluster
    for j = 1:numSample
        contour(x, y, fvnew(:,:,i), "LineWidth",3); hold on;
        contour(x, y, f(:,:,j), '--black'); hold on;
    end
end


