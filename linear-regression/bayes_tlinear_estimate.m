function [glm] = bayes_tlinear_estimate (glm,y)
% Bayesian parameter estimation with MVT-priors for Linear Model
% FORMAT [glm] = bayes_tlinear_estimate (glm,y)
%  
% glm       Input fields are
%
%           .X      [N x p] design matrix of independent variables
%           .b0,.c0 define Gamma prior over noise precision
%           .w0,.B0 define Normal over regression coeffs where covariance
%           is scaled by noise precision. This implies MVT over coeffs.
%
% y         [N x 1] vector containing dependent variable
%
% glm       Output fields are
%
%           .F      Log Model Evidence
%           .bN,.cN,.wN,.BN define posterior MVT
%
% Priors and Posteriors are Multivariate T (MVT) distributions as in [1,2]
%
% [1] Bernardo and Smith, Bayesian Theory, Wiley, 2000. See Appendix A.
%
% [2] W. Penny. Bayesian General Linear Models with T-Priors. 
% Technical Report, Wellcome Trust Centre for Neuroimaging, 2013. 

[N,K]=size(glm.X);
glm.N=N;
glm.K=K;

if ~isfield(glm,'c0') | isempty(glm.c0), glm.c0=1.1; end
if ~isfield(glm,'w0') | isempty(glm.w0), glm.w0=zeros(K,1); end
if ~isfield(glm,'B0') | isempty(glm.B0), glm.B0=eye(K); end
if ~isfield(glm,'b0') | isempty(glm.b0),
    my = mean(y);
    vy = std(y)^2;
    y = y-my;
    glm.b0 = 1/(glm.c0*vy);
end

glm.lambda0=glm.b0*glm.c0;

y=y(:);
NN=length(y);
if ~(NN==glm.N)
    disp('Error in bayes_tlinear_estimate.m: dimensions of X and y don''t match');
    return
end

% Compute Posterior
glm.BN=glm.B0+glm.X'*glm.X;
glm.wN=inv(glm.BN)*(glm.B0*glm.w0+glm.X'*y);

e=y-glm.X*glm.wN;
z=glm.wN-glm.w0;
sse=(1/glm.b0) + 0.5*e'*e + 0.5*z'*glm.B0*z;
glm.bN=1/sse;
glm.cN=glm.c0+0.5*glm.N;
glm.lambdaN=glm.bN*glm.cN;

% Compute Model Evidence
F = gammaln(glm.cN)-gammaln(glm.c0)+glm.c0*log(2/glm.b0)-glm.cN*log(2/glm.bN);
glm.F = F+0.5*(log(det(glm.B0))-log(det(glm.BN)))-0.5*N*log(pi);
