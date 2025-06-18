
clear all
close all

% This script demonstrates MLM using different shrinkage priors for spm_mlm_bayes

% Set dimensionality
N=100;
d=10;
p=5;

% Generate data
x=randn(N,p);
W=randn(p,d);
W(4:5,:)=0; % In this example inputs 4 and 5 are irrelevant
e=2*randn(N,d); % Noise

% Generate outputs
y=x*W+e;

% Global shrinkage prior
options.pr = 'global';
options.verbose=1;
evalc('mlm_global = spm_mlm_bayes (y,x,options);');

% Input-specific shrinkage prior
options.pr = 'input';
evalc('mlm_in = spm_mlm_bayes (y,x,options);');

% Output-specific shrinkage prior
options.pr = 'output'
evalc('mlm_out = spm_mlm_bayes (y,x,options);');

% Compare model evidences
disp('Comparison of shrinkage priors');
disp(sprintf('Log evidence for IS versus global = %1.2f',mlm_in.fm-mlm_global.fm));
disp(sprintf('Log evidence for OS versus global = %1.2f',mlm_out.fm-mlm_global.fm));
disp(sprintf('Log evidence for IS versus OS = %1.2f',mlm_in.fm-mlm_out.fm));

% Choose model for visualisation
mlm=mlm_in;

% Display posterior mean of regression weights - Bayesian
figure
imagesc(mlm.wmean);
colormap gray
colorbar
ylabel('Inputs');
xlabel('Outputs');
title('Bayes Regression Coefficients');

% Display likelihood regression weights (non-regularised)
figure
imagesc(mlm.wml);
colormap gray
colorbar
ylabel('Inputs');
xlabel('Outputs');
title('ML Regression Coefficients');