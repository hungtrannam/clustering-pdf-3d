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

fv = f(:,:,randperm(numSample, numCluster));

fv_copy = fv;

while iter < max_iter
    iter = iter + 1;

    % Compute Wf
    Wf = zeros(numCluster, numSample);
    for j = 1:numSample
        for i = 1:numCluster
            diff = min(fv(:,:,i), f(:,:,j));
            Wf(i, j) = -log(Integration(h, diff, 2)) + 10^(-10);
        end
    end

    [~, IDX] = min(Wf);

    for i = 1:numCluster
        if sum(IDX == i) > 0
            fvnew(:, :, i) = mean(f(:, :, IDX == i), 3);
        end
    end

    ObjFun = sum(min(Wf));
    fprintf('Iteration count = %d, obj. kmeans = %f\n', iter, ObjFun);


    % Check convergence
    
    if  sum(abs(fv - fvnew)) < epsilon
        break;
    end

    % Update for the next iteration
    fv = fvnew;


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


figure(5);
for i = 1:numCluster
    for j = 1:numSample
        contour(x, y, fvnew(:,:,i), "LineWidth",3); hold on;
        contour(x, y, f(:,:,j), '--black'); hold on;
    end
end


