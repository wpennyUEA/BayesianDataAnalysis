
clear all
close all

load pima-dataset

N = size(X,1);
%N = 200;
ind = randperm(N);

K=8;
u = X(ind,1:K);
y = X(ind,9);
[N,K] = size(u);

M.verbose = 1;

% Initial latent space dimension
M.p = 4; 

%M.opt = 'FixedStepSize';

algorithm='BR';
M.lambda=1; % Regularisation parameter
switch algorithm
        
    case 'BR',
        disp('Bayesian Regularization Algorithm');
        M.prune_its_min=1000;
        M = nca_prune (u, y, M);
        
    case 'BP',
        disp('Bayesian Pruning Algorithm');
        M = nca_prune (u, y, M);
        
    otherwise
        disp('Error in demo_nca_synth.m: unknown algorithm');
        return
        
end