function [logbf] = savage_dickey_vector (prior,post,C,h0)
% Savage-Dickey test over vectors for Normal densities
% FORMAT [logbf] = savage_dickey_vector (prior,post,C,h0)
%
% prior.m       mean
% prior.C       covariance
% post.m        posterior mean
% post.C        posterior covariance
% C             contrast matrix of interest
% h0            null hypothesis (default is zero vector)
%
% logbf         Log Bayes Factor in favour of alternative hypothesis

post_mean=C'*post.m;
post_cov=C'*post.C*C;

prior_mean=C'*prior.m;
prior_cov=C'*prior.C*C;

% Enforce exact symmetry - to avoid potential numerical issues
prior_cov=(prior_cov+prior_cov')/2;
post_cov=(post_cov+post_cov')/2;

% Hypothesized value
k=length(prior_mean);
if nargin < 4 | isempty(h0)
    h0=zeros(k,1);
end
    
% For now just use spm_mvNpdf. This can be made more efficient !
%prior_prob=spm_mvNpdf(h0,prior_mean,prior_cov);
%post_prob=spm_mvNpdf(h0,post_mean,post_cov);
%logbf=log(prior_prob+eps)-log(post_prob+eps);

log_prior_prob = logmvNpdf_robust(h0,prior_mean,prior_cov);
log_post_prob = logmvNpdf_robust(h0,post_mean,post_cov);
logbf = log_prior_prob - log_post_prob;

if isinf(logbf) | isnan(logbf)
    keyboard
end
