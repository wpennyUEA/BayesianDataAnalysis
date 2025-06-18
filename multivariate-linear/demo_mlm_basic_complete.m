
clear all
close all

% This script demonstrates how to use spm_mlm_bayes for Bayesian
% Multivariate Linear Modelling

% Set dimensions
N=100;  % Sample number
d=2;    % DV number
p=3;    % IV number

% Generate random data
y=randn(N,d);
x=randn(N,p);

% Set options for Bayesian MLM
options.pr = 'input';   % Group coefficients by input variables
options.verbose = 0;    % Show iteration details during fitting

% Run Bayesian MLM
evalc('mlm = spm_mlm_bayes(y, x, options);');

% Display estimated regression coefficients
figure
imagesc(mlm.wmean);
colorbar
ylabel('Inputs');
xlabel('Outputs');
colormap(gray);

% Display model evidence
disp(sprintf('Model evidence = %1.2f', mlm.fm));