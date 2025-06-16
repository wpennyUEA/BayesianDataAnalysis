
clear all
close all

% ------------------------------------------------------
% Create some data - REPLACE THIS WITH YOUR X, X0 and y
N = 100;
X0 = randn(N,2);
X0 = zmuv(X0); % ensure other variables have zero mean
X = [ones(N,1),X0];

%true_beta = [0,0.3,0.2]';
true_beta = [0,0.3,0.2]';  % change the size of the constant term?

e = 0.1*randn(N,1); % noise
y = X*true_beta + e;

% ------------------------------------------------------
% Set up null and alternative models

% Gamma prior on noise precision with mean=b0*c0; var = (b0^2)*c0;
glm.b0 = 10; glm.c0=0.1; 
glm0 = glm;
glm.X = X;
glm0.X = X0;


% ------------------------------------------------------
% Fit model(s)
glm = bayes_tlinear_estimate (glm,y);

% The "Savage-Dickey" test
logbf = bayes_tlinear_test(glm,[1 0 0],0)

% The "fit-two-models" test
glm0 = bayes_tlinear_estimate (glm0,y);

logBF = glm.F-glm0.F

