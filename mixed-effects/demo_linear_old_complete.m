
clear all
close all

% This script compares Canonical Correlation Analysis (CCA) and
% Nultivariate Linear Models (LML) for predicting simple linear mapping
% between two datasets

% Simulation parameters
rho = 0.5;  % Relationship between CCA and MLM complexity
d1 = 50;    % Dimensionality of dataset 1
d2 = 50;    % Dimensionality of dataset 2

% Number of latent components in CCA
k = ceil(2*d1*d2/(d1+d2)*rho); 
m = 1;  % Number of clusters (1 = standard CCA)

N = 1000;   % Sample number
sigma = 0.1;    % Noise SD
e1 = sigma*randn(d1,N); % Add noise to x1

% Generate synthetic data. x1 is a linear transformation of x2 + noise
x2 = randn(d2,N);
T = randn(d1,d2);   % True transformation matrix
x1 = T*x2 + e1;
% x1 = [x1;ones(1,N)];
% x2 = [x2;ones(1,N)];

% Fit CCA model
options.maxIter = 512;
options.tol = 10^(-6);
tic; cca = vbcca(x1,x2,k,m,options); els=toc    % Variational bayesian CCA

% Predict x1 from x2
con.Gamma1 = eye(d1);
con.Gamma2 = eye(d2);
[x1pred,gamma,p1,T12] = vbcca_cond_subspace (cca,x2,con);
x1_r2 = compute_r2 (x1',x1pred');

% Predict x2 from x1
x2pred = vbcca_cond_context (cca,x1,con);
x2_r2 = compute_r2 (x2',x2pred');

% Plot true vs estimated transformation matrix
% for i=1:d1, subplot(d1,1,i); plot(x1(i,:),'b'); hold on; plot(x1pred(i,:),'r'); end
figure
subplot(2,2,1);
imagesc(T); colorbar
title('True')
subplot(2,2,2);
imagesc(T12{1}); colorbar
title('CCA Estimated')
subplot(2,2,3);
imagesc(T-T12{1}); colorbar
title('Difference')

% Fit MLM model
options.pr='global';
options.ml_only = 0;
options.verbose = 1;
options.pseudo = 1;

% Estimate transformations
tic; txt1 = evalc('mlm_r1 = spm_mlm_bayes(x1'',x2'',options);'); toc % predict region1 data
tic; txt2 = evalc('mlm_r2 = spm_mlm_bayes(x2'',x1'',options);'); toc % predict region2 data

% Predict using learned weights
x1hat = mlm_r1.wmean'*x2;
x2hat = mlm_r2.wmean'*x1;
x1_r2_mlm = compute_r2 (x1',x1hat');
x2_r2_mlm = compute_r2 (x2',x2hat');

% Plot true vs MLM-estimated weights
figure
subplot(2,2,1);
imagesc(T); colorbar
title('True')
subplot(2,2,2);
imagesc(mlm_r1.wmean'); colorbar
title('MLM Estimated')
subplot(2,2,3);
imagesc(T-mlm_r1.wmean'); colorbar
title('Difference')

% Compare prediction R^2 from CCA and MLM
figure
plot(x1_r2);
hold on
plot(x1_r2_mlm,'r');
legend({'CCA','MLM'})
xlabel('x1 variable')
ylabel('R^2')
grid on
figure
plot(x2_r2);
hold on
plot(x2_r2_mlm,'r');
legend({'CCA','MLM'})
xlabel('x2 variable')
ylabel('R^2')
grid on

% Compare model complexity (parameter number)
p_CCA = m*k*(d1+d2);
p_MLM = 2*d1*d2;

disp(sprintf('Number of CCA params = %d',p_CCA));
disp(sprintf('Number of MLM params = %d',p_MLM));


