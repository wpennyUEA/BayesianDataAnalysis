
clear all
close all

% This script uses MCCA models with different numbers of clusters to estimate factor matrices from two data sources
% with a shared latent variable

% Data points
N = 100;

% Generate latent variable
z = randn(1,N);

% Define factor matrices for first two clusters
W1{1} = [0.5 -0.3]'; % for first
W2{1} = [1 2]';
W1{2} = [-0.3 0.5]'; % and second cluster
W2{2} = [2 1]';
% Means 
mu1{1} = [-3 3]'; % for first  
mu2{1} = [3 -3]';
mu1{2} = [2 -2]'; % and second cluster
mu2{2} = [2 2]';

% Observation noise SD
sigma = 0.2;

% Generate data sources - mixing latent variable, cluster factors, means and noise
X1 = []; X2 = [];
for m = 1:2,
    X1new = W1{m}*z + mu1{m}*ones(1,N) + sigma*randn(2,N);
    X2new = W2{m}*z + mu2{m}*ones(1,N) + sigma*randn(2,N);
    X1 = [X1, X1new];
    X2 = [X2, X2new];
end

% Evaluate model evidence for different number of clusters (1-5)
for m = 1:5,
    cca = vbcca (X1,X2,1,m);
    F(m) = cca.F;   % Store model evidence
end

% PLot model evidence vs number of clusters
figure
plot(F);
xlabel('Number of Clusters');
ylabel('Model Evidence');
grid on
title('True number of clusters = 2');

