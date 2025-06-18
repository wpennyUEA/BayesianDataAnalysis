
clear all
close all

% This script applies Bayesian Logistic Regression to the Ionosphere dataset

% Load data
load iono_data.mat

% Labels
t = iono_data(:,35) - 1;
% Select features
x = iono_data(:,[1,3:34]);  % 2nd variable has zero SD, so ignore
[N,P] = size(x);    % Number of samples and features
opt.verbose=1;

% Bayesian Logistic Regression
M = blr_fit (x,t,opt);

% Plot z-scores
figure
plot(M.z(1:end-1));
xlabel('Parameter');
ylabel('Z-score');
grid on
hold on
% Plot significance threshold (95% confidence)
plot([1 P],[1.96,1.96],'r-');
plot([1 P],[-1.96,-1.96],'r-');

% Compare models
model(1).x=[];  % No features, only intercept
model(2).x=x;   % All features
model(3).x=x(:,26); % Only feature 26
% Features with significant z-scores
ind=find(abs(M.z(1:P))>1.96);
model(4).x=x(:,ind);    % Only significant features
name={'const','All','feat26','All-Sig'};

% Bayesian Logistic Regression comparison
F = blr_compare(model,name,t,opt);
