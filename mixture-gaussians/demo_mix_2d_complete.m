
clear all
close all

% This script demonstrates Variational Bayesian GMM in 2D data with a fixed
% number of components

% Load data
load yrep

% Plot raw data points
figure
plot(y(:,1),y(:,2),'x');
hold on

% Assumed number of mixture components (clusters)
m_model=5;
disp('Two-dimensional data with three clusters');
disp(sprintf('Assumed model has %d clusters',m_model));
disp('VB GMM code');

% Variational Bayesian GMM
vbmix=spm_mix(y,m_model);

% Plot estimated means of Gaussian components
for i=1:m_model,
   plot(vbmix.state(i).m(1),vbmix.state(i).m(2),'rx');
end
hold on

% Plot Gaussian mixture components as ellipses
spm_mix_plot2d(vbmix,[-2 12 -2 12],1,'r',0.4,0.5);
set(gca,'FontSize',18);
