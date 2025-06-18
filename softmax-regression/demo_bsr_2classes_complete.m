
clear all
close all

% This script uses Bayesian Softman Regression for classifying two classes
% of 2D data points

% Data points
N=120; 
opt.verbose=1;

% Generate data from mixture model
mix.m=2;    %  Class number
% Means of each Gaussian component
mix.state(1).m=[1,1];
mix.state(2).m=[3,3];
% mix.state(1).m=[1,3];
% mix.state(2).m=[3,1];

% Covariance matrices
for i=1:2,
    mix.state(i).C=eye(2);
    mix.state(i).prior=1/2;
end

% Sample N points from mixture model
[x,label] = spm_samp_mix (mix,N);

% Plot data points
if opt.verbose
    figure
    col={'rx','bx','kx'};
    for i=1:2,
        ind=find(label==i);
        plot(x(ind,1),x(ind,2),col{i});
        hold on
    end
end

% Randomly permute data - avoids order effects
ind = randperm(N);
x = [x(ind,:),ones(N,1)];   % Add bias term
label = label(ind,:);
opt.diag=0; 
opt.verbose=1;

% Bayesian Softmax Regression model
tic; bsr = bsr_fit (x,label,opt); toc

% Plot decision boundary
bsr_plot_boundary (bsr,x);
