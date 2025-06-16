function [bf10] = linregbf(R2,N,p,s)
% Compute Bayes Factor of Linear Regression Model
% FORMAT [bf10] = linregbf(R2,N,p,s)
%
% R2        proportion of variance explained
% N         Number of data points used to fit model
% p         Number of regressors
% s         Expected effect size (default=1)
%
% bf10      Bayes factor in favour of alternative
%
% Null model is regression model with constant term only.
% BF10 -> inf as R2->1 (a desirable consistency property [1])
% 
% [1] Rouder and Morey, Default Bayes Factors for Model Selection 
% in Regression. Multivariate Behavioural Research, 2013.

% More precisely s is the scale parameter of Cauchy Prior on regression
% coefficient (see page 888, [1])
if nargin < 4 | isempty(s)
    s=1;
end

bf10 = integral(@(g) bf_integrand (g,R2,N,p,s),0,Inf);

% -------------------------------------------
function [y] = bf_integrand (g,R2,N,p,s)

% Function to be integrated - equation (11) in [1]
% But take log of each term and exponentiate at end, for numerical stability (as in corrbf.m)

t1 = ((N-p-1)/2)*log(1+g);
t2 = (-(N-1)/2)*log(1+g*(1-R2));
t3 = (-3/2)*log(g);
t4 = log((s*sqrt(N/2))/gamma(1/2));
t5 = (-N*s^2)./(2*g);

y = exp(t1+t2+t3+t4+t5);
