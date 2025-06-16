function [glm] = bayes_glinear_estimate (glm,y)
% Bayesian parameter estimation with Gaussian priors for Linear Model
% FORMAT [glm] = bayes_glinear_estimate (glm,y)
%  
% glm   Input fields are:
%
%       .X      [N x p] design matrix of independent variables
%       .pE     Prior Mean (default is zero vector)
%       .pC     Prior Covariance (default is IID unit variance)
%       .Ce     Error Covariance (default is IID estimated variance)
%
% y     [N x 1] vector containing dependent variable
%
% glm   Output fields are
%
%       .Ep     Posterior Mean
%       .Cp     Posterior Covariance
%       .Rp     Posterior Correlation
%       .z      Z-Scores for each regressor
%       .R2     Proportion of variance explained
%       .pval   Two-sided p-value
%       .F      Log Model Evidence
%
% The above means, covariances and correlations pertain to the
% regression coefficients. 
%
% The prior over regression coefficients is assumed Gaussian.
% Observation Noise Cov, Ce, is assumed known or uncertainty in its estimation
% is not accommodated. The posterior is therefore Gaussian and inferences 
% will be overconfident for small N.

X = glm.X;
y = y(:);

scale=0;
if scale
    X = zmuv(X);
    y = zmuv(y);
    N = size(X,1);
    % Add a column of 1's to the design matrix.
    X = [X,ones(N,1)];
end
[N,p]=size(X);

if ~isfield(glm,'pE') | isempty(glm.pE), pE=zeros(p,1); else pE=glm.pE; end
if ~isfield(glm,'pC') | isempty(glm.pC), pC=eye(p); else pC=glm.pC; end

if ~isfield(glm,'Ce') | isempty(glm.Ce)
    % Compute error variance
    beta=pinv(X)*y;
    err=(y-X*beta);
    s=std(err);
    s = s * (N-1)/(N-p); % Inflate estimate appropriately
    var=s^2;
    glm.Ce=var*eye(N); % IID Error matrix
end
Ce = glm.Ce;

ipC=inv(pC);
iCe=inv(Ce);

% Estimation equations from e.g. Chapter 3 of C Bishop. 
% Pattern Recognition and Machine Learning, Springer, 2006.
iCp=X'*iCe*X+ipC;   % Posterior Precision
Cp=inv(iCp);        
Ep=Cp*(X'*iCe*y+ipC*pE); % Posterior Mean

sp = sqrt(diag(Cp));
z = Ep./sp;
pval = 2*(1-spm_Ncdf(abs(z),0,1));

Rp=Cp./(sp*sp'); 

% Log evidence
yhat = X*Ep;
ey = y-yhat;
ew = pE-Ep;

F = -0.5*spm_logdet(Ce) - 0.5*N*log(2*pi);
F = F - 0.5*trace(ey'*iCe*ey);
F = F - 0.5*ew'*ipC*ew - 0.5*spm_logdet(pC) + 0.5*spm_logdet(Cp);

% Variance explained
ssy = sum(sum(y.^2));
sse = sum(sum(ey.^2));
R2 = (ssy-sse)/ssy;

glm.pE=pE;
glm.pC=pC;
glm.Ep=Ep;
glm.Cp=Cp;
glm.Rp=Rp;
glm.z=z;
glm.pval=pval;
glm.R2=R2;
glm.F=F;
