
clear all
close all

% This script uses a Bayesian Softmax Regression model to classify data
% from a three-class 2D Gaussian mixture model

N=120; % Data points
Dnoise=0; % Number of additional noisy inputs (distractors)
disp(sprintf('Number of spurious predictors added = %d',Dnoise));
opt.verbose=1;

% Mixture model parameters
mix.m=3;
% Means of Gaussian components
mix.state(1).m=[1,1]';
mix.state(2).m=[1,5]';
mix.state(3).m=[3,3]';

% Covariance matrices and equal priors
for i=1:3,
    mix.state(i).C=eye(2);
    mix.state(i).prior=1/3;
end

% Sample N data points
[x,label] = spm_samp_mix (mix,N);

% Plot data points
if opt.verbose
    figure
    col={'rx','bx','kx'};
    for i=1:3,
        ind=find(label==i);
        plot(x(ind,1),x(ind,2),col{i});
        hold on
    end
end

% Randomly permute data - avoids order effects
ind = randperm(N);
x = x(ind,:);
label = label(ind,:);

% Add spurious  predictors (random Gaussian noise)
x = [x,1.65*randn(N,Dnoise),ones(N,1)];

opt.diag=0; % Full covariance model
opt.verbose=1;

% Bayesian Softmax Regression model
tic; bsr = bsr_fit (x,label,opt); toc

% Calculate predictived class probabilities and auxiliary quantities
[y,a] = bsr_output (bsr,x);

% Log Baues Factors and feature relevance indicators
[logbf,z] = bsr_savage_dickey(bsr); % Savage-Dickey
disp(' ');
disp('Log BF, row is class(k), column is feature (d):');
disp(logbf)
disp(' ');
disp('z, row is class(k), column is feature (d):');
disp(z)
disp(' ');
disp('Sum LogBF over classes:');
disp(sum(logbf,1)); % Overall feature importance

% Calculate updated predictions, auxiliary quantities and modified outputs
[y,a,ymod] = bsr_output (bsr,x);

% Display predicted class assignments and decision boundaries
if opt.verbose
    figure
    col={'rx','bx','kx'};
    for i=1:3,
        [tmp,assign]=max(y');   % Maximum predicted probability
        ind=find(assign'==i);
        plot(x(ind,1),x(ind,2),col{i});
        hold on
    end
    % Plot decision boundary
    if Dnoise==0
        % If Dnoise>0 there will be > 2 inputs
        bsr_plot_boundary (bsr,x);
    end
       
end