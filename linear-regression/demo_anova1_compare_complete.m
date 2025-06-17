clear all
close all

% This script simulates data from two groups. It runs a Bayesian One-Way
% ANOVA and compares the Bayes Factor results from 3 approaches - Default,
% Gaussian and T-Distribution

% Name your groups here
group(1).name='Frogs';
group(2).name='Bears';

% Define true means
mu=[3,5];
%mu=[3,3];
mu=mu+100;

% Sample sizes
N=[12,13];

% Name Bayes Factor approaches
approach(1).name='Default';
approach(2).name='Gaussian';
approach(3).name='T';

% Number of simulation repetitions
Reps=20;
for r=1:Reps,
    % Generate simulated data for each group with normal noise around the
    % mean
    group(1).x = randn(N(1),1)+mu(1);
    group(2).x = randn(N(2),1)+mu(2);
    
    % Default Bayes Factor from Rouder 2009
    % implemented in matlab by Schwarzkopf
    [p,stats] = glm_anova1 (group);
    t = sqrt(stats.F);
    BF = t2smpbf(t,N(1),N(2));
    approach(1).logBF(r) = log(BF);
 
    % Using Gaussian method
    approach(2).logBF(r) = bayes_glm_anova1 (group,'G');

    % Using T-Distribution method
    approach(3).logBF(r) = bayes_glm_anova1 (group,'T');

end

% Plot pairwise comparisons of Log Bayes Factors for different approaches
figure
k=1;
for i=1:3,
    for j=i+1:3,
        subplot(2,2,k);
        plot(approach(i).logBF,approach(j).logBF,'x');
        xlabel(approach(i).name);
        ylabel(approach(j).name);
        grid on
        [tmp,ind]=sort(approach(i).logBF);
        hold on
        plot(tmp,tmp,'r-');
        k=k+1;
    end
end
