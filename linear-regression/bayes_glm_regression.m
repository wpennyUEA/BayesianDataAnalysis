function [logBF,glm] = bayes_glm_regression (X,y,method)
% Bayesian GLM regression
% FORMAT [logBF,glm] = bayes_glm_regression (X,y,method)
%
% X                 [N x K] Design matrix (Independent Variables).
%                   Kth column must be vector of 1's.
% y                 [N x 1] Dependent Variable
% method            'G'     Multivariate Gaussian priors and posteriors
%                   'T'     Multivariate T (MVT) priors and posteriors
%
%                   Log Bayes Factors in favour of alternative hypothesis:
% logBF(1)          Fit Full and Reduced Model separately (Multi-Fitting)
% logBF(2)          Savage-Dickey
% logBF(3)          Default Bayes Factors
%
% glm               Linear Model parameters
%
% Single variance component assumed

[N,K] = size(X);

% Check Kth column is vector of 1's
d=X(:,K)-ones(N,1);
if sum(d.^2) > 0
    disp('Error in bayes_glm_regression: last column of design matrix must be column of 1''s.');
    return
end
    
% Centre design matrix - as in Rouder and Morey (2013)
for p=1:K-1,
    X(:,p)=X(:,p)-mean(X(:,p));
end

% Subtract data mean
my = mean(y);
vy = std(y)^2;
y = y-my;

% Prior precision/variance matrices
XK = X(:,1:K-1);
B0 = (XK'*XK)/N;
B0 = blkdiag(B0,1);

% Rouder and Morey (2013) set the prior covariance of regression
% coefficients to be proportional to inv(B0). Hence choice of the
% above term.

switch method
    case 'T',
        % Single variance component assumed - intrinsic to MVT prior
        
        % c0 must be larger than 1 for finite variance
        % 2*c0 is the degrees of freedom of T-prior
        % (c0=1/2 for Cauchy)
        glm.c0=1.1;
        glm.B0=B0; 
        glm.w0=zeros(K,1);
        
        % Scale observation noise level to variance of dependent variable.
        % This then changes scale of prior variance over regression
        % coefficients through T-prior.
        glm.b0 = 1/(glm.c0*vy);
        
        % Define null model
        glm_null = glm;
        glm_null.w0 = glm.w0(K);
        glm_null.B0 = glm.B0(K,K);

    case 'G',
        glm.pE=zeros(K,1);
        
        % Single variance component assumed
        glm.Ce=vy*eye(N); % IID Error matrix
        
        % Scale prior variance by data variability
        C0 = inv(B0);
        glm.pC=vy*C0;
        
        % Define null model
        glm_null = glm;
        glm_null.pE = glm.pE(K);
        glm_null.pC = glm.pC(K,K);
        
    otherwise
        disp('Unknown method in bayes_glm_regression.m');
        return
end

glm.X = X;
glm_null.X = glm.X(:,K);

glm = bayes_linear_estimate (glm,y,method);

% Savage-Dickey
c = eye(K);
logBF(2) = bayes_linear_test (glm,c,method);

% Fit null model
glm_null = bayes_linear_estimate (glm_null,y,method);
logBF(1) = glm.F - glm_null.F;

% Default Bayes Factor
R2 = glm.R2;
bf10 = linregbf(R2,N,K-1);
logBF(3) = log(bf10);

