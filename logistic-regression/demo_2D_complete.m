
clear all
close all

% This script demonstrate Bayesian Logistic Regression to classify
% two-class Gaussian data in 2D. 

N=120; % Data points
opt.verbose=1;

mix.m=2;    % Class number
% Mean vectors
mix.state(1).m=[1,1];
mix.state(2).m=[3,3];
% Covariance matrices and equal class priors
for i=1:2,
    mix.state(i).C=eye(2);
    mix.state(i).prior=1/2;
end

% Generate N samples
[x,label] = spm_samp_mix (mix,N);

% Plot data
if opt.verbose
    figure
    col={'rx','bx'};
    for i=1:2,
        ind=find(label==i);
        plot(x(ind,1),x(ind,2),col{i});
        hold on
    end
end

% Randomly permute data - avoids order effects
ind = randperm(N);
x = x(ind,:);
label = label(ind,:);

opt.diag=0; % Full covariance
opt.verbose=1; 

% Labels for logistic regression
t = label-1;

% Bayesian Logistic Regression
M = blr_fit (x,t,opt);

% Plot decision boundary
blr_plot_boundary (M,x);
