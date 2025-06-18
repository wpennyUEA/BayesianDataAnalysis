
clear all
close all

% This script uses Bayesian CCA to estimate the factor matrices of two
% independent 2D data sources

% Data points
N = 500;

% Generate independent data sources
X1 = randn(2,N);
X2 = randn(2,N);

% Options for CCA algorithm
options.maxIter = 512;
options.tol = 10^(-5);

% Variational Bayesian CCA
cca = vbcca (X1,X2,1,1,options);    % One latent variable

% Display estimated factor matrices and means for both data sources
disp('Data source 1:');
disp('Estimated W');
disp([cca.W1{1}]);
disp('Estimated mean');
disp([cca.mu1{1}]);
disp('Data source 2:');
disp('Estimated W');
disp([cca.W2{1}]);
disp('Estimated mean');
disp([cca.mu2{1}]);

% Plot model evidence over iterations
figure
plot(cca.Fhist);
ylabel('Model Evidence');
xlabel('Number of iterations');
grid on

% VB CCA  with null model (no latent variable)
options.null = 1;
cca_null = vbcca (X1,X2,1,1,options);

% Calculate log Bayes Factor comparing model and null
logBF_alt = cca.F - cca_null.F

