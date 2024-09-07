clear; close all;
m = 2;  
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
U0 = rand(numCluster, numSample);
U0 = U0 ./ sum(U0);
fv = (f * U0.^m') ./ sum(U0.^m, 2)';
fv_copy = fv;

while iter < max_iter
    iter = iter + 1;

    for j = 1:numSample
        for i = 1:numCluster
            Wf(i, j) = Integration(h, abs(fv(:, i) - f(:, j)), 1) + 10^(-10);
        end
    end

    for i = 1:numCluster
        for j = 1:numSample
            numerator = 1 / (Wf(i, j)) ^ (2/(m - 1));
            denominator = 0;
            for k = 1:numCluster
                denominator = denominator + 1 / (Wf(k, j)) ^ (2/(m - 1));
            end
            Unew(i, j) = numerator / denominator;
        end
    end

    fvnew = (f * Unew.^m') ./ sum(Unew.^m, 2)';


    ObjFun = sum(sum((Unew) .* Wf));
    fprintf('Iteration count = %d, obj. ifcm = %f\n', iter, ObjFun);

    if norm(fv - fvnew) < epsilon
        break;
    end

    fv = fvnew;
    U = Unew;

    [~,IDX] = max(Unew);
    for i = 1:numCluster
        plot(x, f(:,IDX == i),'Color', color{i}, 'LineStyle',':', 'LineWidth', 2.5, 'DisplayName', 'Density Function Vectors'); hold on;
        % plot(x, fv_copy(:,i), 'LineWidth', 2, 'Color', color{i});
        plot(x, fv(:,i), 'LineWidth', iter/1, 'Color', color{i}, 'DisplayName', 'Prototypes');
        box off;
        set(gca, 'YColor', 'none');
    end

end

figure;
colormap(abyss);
b = bar3(Unew);

for k = 1:length(b)
    b(k).CData = b(k).ZData;
    b(k).FaceColor = 'interp'; 
end

colorbar;
