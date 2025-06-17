
% This script compares Bayesian/variational PCA to classical maximum likelihood (ML)
% PCA for dimensionality reduction

% Load image
x=imread('Alan.jpg','jpg');

% Crop
startx=100;
starty=300;
N=256; % Size of square cropped image (pixels)

% Convert green chennel to double precision
xg=double(x(startx:startx+N-1,starty:starty+N-1,2));

% Centre image data
xg=xg-mean(mean(xg));

% Set maximum latent space dimensionality (number of components)
q=100;

% Run Variational PCA
pca=spm_vpca(xg,q);

% Posterior mean of weights: Bayesian estimate
figure; imagesc(pca.M_w); colormap gray; title('Bayes estimate');
% ML weights: Maximum Likelihood estimate
figure; imagesc(pca.ml.W(:,1:q)); colormap gray; title('ML estimate');

% Plot negative Free Energy over iterations
figure
plot(pca.Fm_evol);
xlabel('Iterations');
ylabel('Neg. Free Energy');

% Plot eigenspectrum
figure
plot(pca.ml.lambda);
title('Eigenspectrum');

% Plot the learned prior precisionfor each latent factor
figure
plot(pca.mean_alpha);
title('Prior precision of factors');
