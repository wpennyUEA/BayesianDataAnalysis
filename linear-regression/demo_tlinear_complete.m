
clear all
close all

% This script compares the Savage-Dickey and Model comparison Bayes factor
% approaches to Bayesian hypthesis testing

% Add data here
N = 100;    % Sample size
X0 = randn(N,2);    % Two regressors
X0 = zmuv(X0); % Zero mean, unit-variance normalisation
X = [ones(N,1),X0]; % Full design matrix with intercept

% Define true beta weights (intercept and coefficients)
%true_beta = [0,0.3,0.2]';
true_beta = [0,0.3,0.2]';  % Can change the intercept (0 here)

% Generate response variable with Gaussian noise
e = 0.1*randn(N,1); 
y = X*true_beta + e;

% Set up null and alternative models
% Prior for the noise precision (inverse variance) with mean=b0*c0; var = (b0^2)*c0;
glm.b0 = 10; glm.c0=0.1; 
glm0 = glm;  % Null model - excludes some predictors
glm.X = X;   % Alternative model - includes all predictors
glm0.X = X0; % Null model - excludes intercept

% Fit alternative model (includes intercept and X0)
glm = bayes_tlinear_estimate (glm,y);

% Test intercept term using Savage-Dickey method
logbf_savage_dickey = bayes_tlinear_test(glm,[1 0 0],0)

% Fit null model (without intercept)
glm0 = bayes_tlinear_estimate (glm0,y);

% Compare models by subtracting log evidences
logbf_two_model = glm.F-glm0.F

