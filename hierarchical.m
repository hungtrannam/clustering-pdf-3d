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


for j = 1:numSample
    for i = 1:numSample
        Wf(i, j) = Integration(h, abs(f(:, i) - f(:, j)), 1) + 10^(-10);
    end
end


Z = linkage(Wf, "ward");

[H, T, outperm] = dendrogram(Z ,'ColorThreshold','default');
