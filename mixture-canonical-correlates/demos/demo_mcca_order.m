
clear all
close all

disp('Create two 2D data sources with single latent variable and two clusters');
disp('Estimate factor matrices using MCCA');

% Number of data points
N = 100;

% Latent variable
z = randn(1,N);

% Factor matrices  
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

% Generate data sources
X1 = []; X2 = [];
for m = 1:2,
    X1new = W1{m}*z + mu1{m}*ones(1,N) + sigma*randn(2,N);
    X2new = W2{m}*z + mu2{m}*ones(1,N) + sigma*randn(2,N);
    X1 = [X1, X1new];
    X2 = [X2, X2new];
end

% Call MCCA algorithm
for m = 1:5,
    cca = vbcca (X1,X2,1,m);
    F(m) = cca.F;
end

figure
plot(F);
xlabel('Number of Clusters');
ylabel('Model Evidence');
grid on
title('True number of clusters = 2');

