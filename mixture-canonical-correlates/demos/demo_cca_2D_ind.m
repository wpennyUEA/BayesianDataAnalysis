
clear all
close all

disp('Create two independent 2D data sources');
disp('Estimate factor matrices using CCA');

% Number of data points
N = 500;

% Generate data sources
X1 = randn(2,N);
X2 = randn(2,N);

% Call CCA algorithm
options.maxIter = 512;
options.tol = 10^(-5);
cca = vbcca (X1,X2,1,1,options);

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

figure
plot(cca.Fhist);
ylabel('Model Evidence');
xlabel('Number of iterations');
grid on

% Call CCA algorithm with null model option
options.null = 1;
cca_null = vbcca (X1,X2,1,1,options);

logBF_alt = cca.F - cca_null.F

