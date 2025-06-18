
clear all
close all

% This script demonstrates CVA on simulated data with a single latent
% factor. It compares estimated canonical vectors to the true ones

% Dimensionality of both datasets
d1=3;
d2=5;
N=30;

% Max posible number of canonical vectors
%m=min([d1,d2])
m=1;

% Generate true factor loading matrices
W1=10*randn(d1,m);
W2=10*randn(d2,m);

% Observation noise SD
sig=0.01;

% Generate noise
E1=sig*randn(d1,N);
E2=sig*randn(d2,N);

% Generate latent sources and observed data matrices
if m==0
    X1=E1;
    X2=E2;
else
    % Generate latent factor time courses
    Z=randn(m,N);
    % Generate observed data (factor loadings times sources + noise)
    X1=W1*Z+E1;
    X2=W2*Z+E2;
end

% Bayesian CVA model
CVA = spm_cva_prob (X1,X2);

% Absolute values of true and estimated canonical vectors for both datasets
disp('True');
abs(W1)
disp('Estimated');
abs(CVA.W1)
disp('True');
abs(W2)
disp('Estimated');
abs(CVA.W2)
