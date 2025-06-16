
clear all
close all

disp('Create two 2D data sources from CCA model with single latent variable');
disp('and known factor matrices.');
disp('Estimate factor matrices using CCA');

% Number of data points
N = 500;

% Latent variable
z = randn(1,N);

% Factor matrices
W1 = [5 -3]';
W2 = [1 2]';

% Means
mu1 = [1 1]';
mu2 = [3 3]';

% Observation noise SD
sigma = 0.1;

% Generate data sources
X1 = W1*z + mu1*ones(1,N) + 0.01*sigma*randn(2,N);
X2 = W2*z + mu2*ones(1,N) + sigma*randn(2,N);

% Call vbCCA algorithm
options.maxIter = 512;
options.tol = 10^(-5);
cca = vbcca (X1,X2,1,1,options);

% Call CCA algorithm with null model option
options.null = 1;
cca_null = vbcca (X1,X2,1,1,options);

% Call probCCA algorithm
cca_prob = probCCA (X1,X2,1);

disp('Data source 1:');
disp('True W      VB-W    ML-W');
disp([W1,cca.W1{1},cca_prob.W1{1}]);
disp('True mean    VB-mean    ML-mean');
disp([mu1,cca.mu1{1},cca_prob.mu1{1}]);

disp('Data source 2:');
disp('True W       VB-W     ML-W');
disp([W2,cca.W2{1},cca_prob.W2{1}]);
disp('True mean    VB-mean    ML-mean');
disp([mu2,cca.mu2{1},cca_prob.mu2{1}]);

disp('Obs noise covariance for data source 1:')
disp('True:')
disp(sigma^2*eye(2))
disp('VB:')
disp(cca.C1{1})
disp('ML:')
disp(cca_prob.C1{1})

figure
plot(cca.Fhist);
ylabel('Model Evidence');
xlabel('Number of iterations');
grid on

logBF_alt = cca.F - cca_null.F

figure
bar(cca.Fdecomp.term);
set(gca,'XTickLabel',cca.Fdecomp.name);
grid on
ylabel('Energies')
Ne = length(cca.Fdecomp.term);
hold on
plot([0 Ne],cca.F*ones(1,2),'r-');

