
clear all
close all

% This script tests linear hypotheses in a GLM with two predictors, using
% contrast matrices. It also compares the true effect with the estimated
% effect

% Number of observations
N=12;
X=randn(N,2);   % Design matrix
e=randn(N,1);   % Gaussian noise
  
% Choose example case
ex=1;
switch ex,
    case 1,
        disp('Example 1:');
        beta = [0 0]';  % No true effect
        
    case 2,
        disp('Example 2:');
        beta = [3 2]';  % True effects on both regressors
end

% Generate outcome variable using linear model
y = X * beta + e;

% Define two linear hypotheses to test
effect(1).c = [1 0; 0 1];
effect(2).c = [1 -1]';

% Test each effect
for i=1:length(effect)
    stats = glm_test_hypothesis (X,y,effect(i).c);
    % Display results
    disp(' ');
    disp('True effect:');
    disp(effect(i).c'*beta);
    disp(' ');
    disp('Estimated effect:');
    disp(stats.effect);
    disp(sprintf('F=%1.2f, df=(%1.2f,%1.2f), p=%1.6f', stats.F,stats.df(1),stats.df(2),stats.p));
    disp(' ');
end


