function [logbf] = bayes_tlinear_test (glm,c,a)
% Bayesian hypothesis testing for Linear Model with T prior/posterior
% FORMAT [logbf] = bayes_tlinear_test (glm,c,a)
%
% glm       output of bayes_tlinear_estimate.m
% c         [K x P] contrast matrix where P is number of regression
%           coefficients and K is dimension of contrast
% a         value of contrast under null hypothesis 
%
% logbf     Log Bayes Factor in favour of alternative hypothesis
%
% Bayes Factors are computed using the Savage-Dickey approximation
% based on the Multivariate T (MVT) Posterior Distribution

[K,P] = size(c);
Pchk=length(glm.w0);
if ~(P==Pchk), 
    disp('Contrast matrix and glm coefficients have inconsistent dimensions'); 
    return; 
end

prior_mean = c*glm.w0;
var = c*inv(glm.B0)*c';
prior_precision = glm.lambda0*inv(var);
prior_df = glm.c0*2;

post_mean = c*glm.wN;
var = c*inv(glm.BN)*c';
post_precision = glm.lambdaN*inv(var);
post_df = glm.cN*2;


% For now just use spm_mvtpdf. This can be made more efficient !
prior_prob = spm_mvtpdf(a,prior_mean,prior_precision,prior_df);
post_prob = spm_mvtpdf(a,post_mean,post_precision,post_df);

logbf = log(prior_prob)-log(post_prob);
