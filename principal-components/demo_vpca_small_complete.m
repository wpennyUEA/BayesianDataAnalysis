
close all
clear all

% This script demonstrates variational PCA using synthetic data. 
% Data is generated from a known PCA generative model: y_n = W x_n + mu + e_n
% where W is d*q - first example in Bishop's VPCA paper

N=100;  % Observations
d=10;   % Data dimensionality
%q=d-1;
q=3;    % Latent imensionality

% Generate orthogonal latent directions
W=randn(d,q);
W=orth(W);

% Generate latent sources
x=randn(q,N);
sd_x=diag([5,4,3,2,1,1,1,1,1]);% Standard deviations
x=sd_x(1:q,1:q)*x;

% Generate Gaussian sensor noise
e=randn(d,N);
% Generate constant mean offset
mu=ones(d,1)*ones(1,N);

% Final observed data
%e=zeros(d,N);
t=W*x+mu+e;

% Run variational PCA
pca=spm_vpca(t);

% Visualise estimated factor matrices
figure; imagesc(pca.M_w); colormap gray; title('Bayes estimate');
figure; imagesc(pca.ml.W); colormap gray; title('ML estimate');

% Plot convergence of variational free energy
figure
plot(pca.Fm_evol);
xlabel('Iterations');
ylabel('Neg. Free Energy');

% Plot eigenspectrum
figure
plot(pca.ml.lambda);
title('Eigenspectrum');

% Plot inferred prior precision on latent factors
figure
plot(pca.mean_alpha);
title('Prior precision of factors');

% Project data into latent space using the top components
disp(' ');
disp('To recover eg 4 hidden sources')
disp('project data, t, onto first 4 columns of factor matrix:');
disp('W=pca.M_w(:,1:4);xhat=W''*t;');
W=pca.M_w(:,1:4);
xhat=W'*t;  % Reconstructed latent sources

% Estimate Bayesian data covariance matrix
disp(' ');
disp('A Bayesian estimate of the data covariance matrix');
disp('is given by:');
disp('obs_noise_var=(1/pca.mean_tau); CBayes=W*W''+obs_noise_var*eye(pca.d)');
obs_noise_var=(1/pca.mean_tau);
CBayes=W*W'+obs_noise_var*eye(pca.d);

% Compare to the empirical covariance of the data
C=cov(t');

% Plot covariances and comparison
figure
subplot(2,2,1);
imagesc(C);
colormap gray
colorbar
title('Data Covariance using Cov');
subplot(2,2,3);
imagesc(CBayes);
colormap gray
colorbar
title('Data Covariance using Bayes VPCA');
subplot(2,2,2);
imagesc(C-CBayes);
colormap gray
colorbar
title('Difference');

