clear all
close all

% This is comparing the weights of the chicks that are fed with either
% 'Horsebean' or 'Linseed' using Bayesian One-Way ANOVA, with two methods -
% Gaussian and T-Distribution

% Name your groups here
group(1).name='Horsebean';
group(2).name='Linseed';

% Name the apporiaches for Bayes Factor calculations
approach(1).name='Default';
approach(2).name='Gaussian';
approach(3).name='T';

% Chick Weights Example Data from "R" 
group(1).x = [179,160,136,227,217,168,108,124,143,140];
group(2).x = [309,229,181,141,260,203,148,169,213,257,244,271];

% Calculate sample sizes
N1 = length(group(1).x);
N2 = length(group(2).x);

% Default Bayes Factor from Rouder
% implemented by Schwarzkopf
% glm_anova1 performs a simple one-way ANOVA
[p,stats] = glm_anova1 (group);
t = sqrt(stats.F);  % Convert F statistic to t statistic
DBF = t2smpbf(t,N1,N2); % Calculate Bayes Factor

% Calculate Log Bayes Factor using two methods - Gaussian and
% T-Distribution methods
logBFG = bayes_glm_anova1 (group,'G');
logBFT = bayes_glm_anova1 (group,'T');

% Display results
disp(' ');
disp('R Bayes Factor Package:');
disp('BF=5.98:');
disp(' ');
disp('Matlab Default Bayes Factor:');
disp(DBF);
disp('Matlab Default Log Bayes Factor:');
log(DBF)

% Using Gaussian method
disp('Gaussian LogBF:');
disp(logBFG);

% Using T-Distribution method
disp('T-Dist LogBF:');
disp(logBFT);




