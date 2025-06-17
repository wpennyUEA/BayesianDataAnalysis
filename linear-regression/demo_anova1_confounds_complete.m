
close all
clear all

% This script performs a Bayesian One-Way ANOVA (GLM approach) on three groups, and
% includes confounds (covariates) in the analysis

% Name your groups here
group(1).name='Frogs';
group(2).name='Bears';
group(3).name='Bats';

% Add your data and sample sizes in here
N1=12; N2=13; N3=9;
group(1).x = randn(N1,1)+3;
group(2).x = randn(N2,1)+4;
group(3).x = randn(N3,1)+2;

% Choose Bayesian method - G for Gaussuan or T for T-distribution
%method ='G';
method = 'T';

disp(' ');
disp('Without confounds:');
% Run GLM ANOVA without confounds
[p,stats] = glm_anova1 (group,[],1);

% Calculate Log Bayes Factor
logBF = bayes_glm_anova1 (group,method,[],1);

% Check with Matlab stats toolbox
%[p_check,stats_check]=my_anova1(group);

% Add confounds here
disp(' ');
disp('With confounds:');
Xc = randn(N1+N2+N3,3);

% Run GLM ANOVA with confounds
[p,stats] = glm_anova1 (group,Xc,1);

% Calculate Log Bayes Facotr
logBF = bayes_glm_anova1 (group, method, Xc, 1);
