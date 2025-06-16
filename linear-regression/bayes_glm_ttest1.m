function [logbf,logbf_jzs,t,p] = bayes_glm_ttest1 (y)
% One sample t-test using Bayesian regression code with T priors/posts
% FORMAT [logbf,logbf_jzs,t,p] = bayes_glm_ttest1 (y)
%
% y         [N x 1] data vector
%
% logbf     Log Bayes Factor in favour of alternative
% logbf_jzs Same but from Default Bayes Factor (JZS) code
% t         t statistic
% p         p-value

y=y(:);
N=length(y);
glm.X=ones(N,1);

% c0 must be larger than 1 for proper density
% 2*c0 is the degrees of freedom of T-prior
glm.c0=1.1; 
glm.B0=1;
glm.w0=0;

% scale observation noise level to data
% this then changes scale of prior variance over regression
% coefficients through T-prior
glm.b0 = 1/(glm.c0*std(y)^2);

glm = bayes_linear_estimate (glm,y,'T');
logbf = bayes_linear_test (glm,1,'T',0);

[h,p,ci,stats] = ttest(y);
t=stats.tstat;
bf10 = t1smpbf(t,N);
logbf_jzs = log(bf10);