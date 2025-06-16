function [logbf] = bayes_glinear_test (glm,c,a)
% Bayesian hypothesis testing for Linear Model with Gaussian priors
% FORMAT [logbf] = bayes_glinear_test (glm,c,a)
%
% glm       output of bayes_linear_estimate.m
% c         [K x P] contrast matrix where P is number of regression
%           coefficients and K is dimension of contrast
% a         value of contrast under null hypothesis 
%
% logbf     Log Bayes Factor in favour of alternative hypothesis
%
% Bayes Factors are computed using the Savage-Dickey approximation
% based on the Multivariate Gaussian Posterior Distribution

[K,P] = size(c);
Pchk=length(glm.pE);
if ~(P==Pchk), 
    disp('Contrast matrix and glm coefficients have inconsistent dimensions'); 
    return; 
end

prior_mean=c*glm.pE;
post_mean=c*glm.Ep;

prior_cov = c*glm.pC*c';
post_cov = c*glm.Cp*c';

% Enforce symmetry - to avoid numerical issues
prior_cov=(prior_cov+prior_cov')/2;
post_cov=(post_cov+post_cov')/2;

% For now just use spm_mvNpdf. This can be made more efficient !
prior_prob=spm_mvNpdf(a,prior_mean,prior_cov);
post_prob=spm_mvNpdf(a,post_mean,post_cov);

logbf=log(prior_prob)-log(post_prob);