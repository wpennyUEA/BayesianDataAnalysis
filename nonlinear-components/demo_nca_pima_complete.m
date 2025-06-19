
clear all
close all

% This script uses Nonlinear Component Analysis (NCA) with Bayesian
% Regularisation to reduce the dimensionality of the 8D input data to 4D latent space

% Load data
load pima-dataset

N = size(X,1);  % Sample size
%N = 200;
ind = randperm(N);  % Random permutation of data indices

K=8;
u = X(ind,1:K);  % Input features
y = X(ind,9);    % Target labels
[N,K] = size(u); % Update dimensions of input

M.verbose = 1;

% Desired latent space dimension
M.p = 4; 

%M.opt = 'FixedStepSize';

% Bayesian Regularisation
algorithm='BR';
M.lambda=1; % Regularisation parameter
switch algorithm
        
    case 'BR',
        disp('Bayesian Regularization Algorithm');
        M.prune_its_min=1000;   % Minimum number of pruning iterations
        M = nca_prune (u, y, M);    % Run NCA
        
    case 'BP',
        disp('Bayesian Pruning Algorithm');
        M = nca_prune (u, y, M);    % Run NCA with pruning
        
    otherwise
        disp('Error in demo_nca_synth.m: unknown algorithm');
        return
        
end