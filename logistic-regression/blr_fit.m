function [M] = blr_fit (X,t,opt)
% Bayesian Logistic Regression
% FORMAT [M] = blr_fit (X,t,opt)
% 
% X     [N x (P-1)] feature vector, (P-1)-dim features, N samples
% t     [N x 1] binary class labels (0/1)
%
% M     .m      [P x 1] posterior mean of parameters
%       .R      [P x P] posterior precision of parameters
%       .z      [P x 1] posterior z-score for parameters
%       .F      model evidence
%
% This routine will add a column of 1's to the design matrix
% making P parameters in total

if nargin < 3 | ~isfield(opt,'verbose'), verbose=0; else verbose=opt.verbose; end
if nargin < 3 | ~isfield(opt,'tol'), tol=1e-4; else tol=opt.tol; end
if nargin < 3 | ~isfield(opt,'alpha_init'), alpha_init=1; else alpha_init=opt.alpha_init; end

Ncheck = length(t);
if isempty(X)
    N = Ncheck; P = 1;
else
    [N,P1] = size(X);
    P = P1 +1;
end
X = [X,ones(N,1)];

if ~(N==Ncheck),
    disp('Error in blr_fit: dimensions of x and t don''t match');
    keyboard
    return
end

M.R0 = eye(P);
M.C0 = inv(M.R0);
M.m0 = zeros(P,1);
Xt = X';

%w = spm_normrnd(M.m0,M.C0,1);
w = zeros(P,1);

maxits = 64;
for it = 1:maxits,
    [J(it),y] = blr_logjoint (w,M,X,t);
    if it > 2,
        if verbose,
            fprintf('Iteration %d, Log Joint = %1.6f \n',it,J(it));
        end
    
        delta_J = (J(it)-J(it-1))/J(it-1);
        if abs(delta_J) < tol, break; end
    end
    
    % Gradient
    g = Xt*(y-t)+M.R0*(w-M.m0);
  
    % Hessian
    V = diag(y.*(1-y));
    H = Xt*V*X+M.R0;
    iH = inv(H);
    
    wold = w;
    dw = - iH*g;
    
    Jlast=J(it);
    Jnew=-Inf;
    alpha=alpha_init;
    for i=1:10,
        w = wold + alpha*dw;
        Jnew = blr_logjoint(w,M,X,t);
        if Jnew < Jlast
            alpha=alpha/2;
        else
            break
        end
    end
    if Jnew < Jlast,
        w=wold;
        break;
    end
    Jlast=Jnew;
    
end

% Posterior Distribution
M.m = w;
M.R = H;
M.C = iH;
M.z = M.m./sqrt(diag(M.C));

% Model evidence
J = blr_logjoint (w,M,X,t);
M.F = J - 0.5*spm_logdet(M.R) + 0.5*P*log(2*pi);