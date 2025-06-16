function [bsr] = bsr_fit (x,label,opt,bsr)
% Fit Bayesian Softmax Regression model
% FORMAT [bsr] = bsr_fit (x,label,opt,bsr)
%
% x     [N x d]      inputs (last column must be 1's)
% label [N x 1]      categories, with label(n) in {1,2,..K}
% opt   .diag        = 1 to use diagonal posterior precision matrix
%                    (and hence no matrix inversion - this may be better for large d)
%       .verbose
%       .winit       initial values of regression coefficients 
%       .r0 prior    prior precision on regression coeffs
%       .alpha_init  initial step size parameter
%       .tol         convergence criterion
%       .hier        use hierarchical prior for new categories
%       .gamma       [N x 1] sample weighting vector (used for mixture models)
%       .S           [Ns x N] selection matrix (used for ignoring task-irrelevant samples)
%
% bsr   .class(k).w0 prior mean provided
%       .class(k).R0 prior precision provided
%       .labels      list of class labels so far encountered
%
% bsr   
%       .class(k).w0  default prior means for class k reg coeffs
%       .class(k).R0  default prior precisions for class k reg coeffs
%       .class(k).iR0 prior cov
%       .class(k).w   posterior means for class k reg coeffs
%       .class(k).R   posterior precisions for class k reg coeffs
%       .class(k).iR  post cov
%       .labels       list of class labels so far encountered
%
% The .gamma and .S fields are used to embed BSR in Cluster-Based Models
% No need to use otherwise

% Default parameters
if nargin < 3 | ~isfield(opt,'diag'), diagpost=0; else diagpost=opt.diag; end
if nargin < 3 | ~isfield(opt,'verbose'), verbose=0; else verbose=opt.verbose; end
if nargin < 3 | ~isfield(opt,'winit'), winit=[]; else winit=opt.winit; end
if nargin < 3 | ~isfield(opt,'r0'), r0=10; else r0=opt.r0; end
if nargin < 3 | ~isfield(opt,'alpha_init'), alpha_init=1; else alpha_init=opt.alpha_init; end
if nargin < 3 | ~isfield(opt,'max_its'), max_its=64; else max_its=opt.max_its; end
if nargin < 3 | ~isfield(opt,'tol'), tol=1e-4; else tol=opt.tol; end
if nargin < 3 | ~isfield(opt,'gamma'), weighting=0; else gamma=opt.gamma; weighting=1; end
if nargin < 3 | ~isfield(opt,'S'), selection=0; else S=opt.S; selection=1; end

[N,d]=size(x);
U=x;
Ut=U';

if sum(abs(U(:,end)-ones(N,1))) > 0
    disp('Error in bsr_fit.m: last column of x must be 1s');
    return
end


if nargin == 4,
    class=bsr.class; 
    Kold=length(class);
    datalabels=unique(label);
    Kdata=length(datalabels);
    newlabels=setdiff(datalabels,bsr.labels);
    labels=[bsr.labels;newlabels(:)];
    K=Kold+length(newlabels);
    
    % Initialise regression coefficients
    for k=1:Kold,
        class(k).w = winit(:,k);
    end
    
    if opt.hier
        hier = hier_fit (bsr);
    end
    
    for k=Kold+1:K,
        disp(sprintf('New category %d :',labels(k)));
        if opt.hier
            disp('Using Hierarchical Prior');
            class(k).w0=hier.mu;
            class(k).R0=hier.Rexp;
        else
            disp('Using Shrinkage Prior');
            class(k).w0=zeros(d,1);
            class(k).R0=r0*eye(d);
        end
        
        iR0 = inv(class(k).R0);
        class(k).iR0 = iR0;
        class(k).w = spm_normrnd(class(k).w0,iR0,1);
    end
else
    labels=unique(label);
    K=length(labels);
    for k=1:K,
        class(k).w0=zeros(d,1);
        class(k).R0=r0*eye(d);
        
        if isempty(winit)
            % Sample from prior
            iR0 = inv(class(k).R0);
            class(k).w = spm_normrnd(class(k).w0,iR0,1);
        else
            class(k).w = winit(:,k);
        end
    end
end

% Create 1-of-K label vectors 
for k=1:K,
    c(:,k)=zeros(N,1);
    ind=find(label==labels(k));
    c(ind,k)=1;
    class(k).iR0 = inv(class(k).R0);
end
    
% To embed BSR in CBM ---------------------------------------------
if selection
    U=S*U;
    Ut=Ut*S';
    c=S*c; 
    if weighting
        gamma=S*gamma;
    end
end
%-------------------------------------------------------------------

for it=1:max_its,
   
    bsr.class=class;
    [y,a] = bsr_output(bsr,U);
    if weighting
        J(it) = bsr_logjoint(bsr,U,c,gamma);
    else
        J(it) = bsr_logjoint(bsr,U,c);
    end
    if verbose,
        fprintf('Iteration %d, Log Joint = %1.6f \n',it,J(it));
    end
    if it > 2,
        delta_J = (J(it)-J(it-1))/J(it-1);
        if abs(delta_J) < tol, break; end
    end
    
    Jlast=J(it);
    for k=1:K,
        
        % Posterior Precision
        Pk=diag(y(:,k).*(1-y(:,k)));
        
        Upre = Ut;
        % To embed BSR in CBM ---------------------------------------------
        if weighting
            Pk = Pk.*diag(gamma);
            Upre = Upre*diag(gamma);
        end
        %------------------------------------------------------------------
        
        class(k).R = Ut*Pk*U+class(k).R0;
        if diagpost
            iR = diag(1./diag(class(k).R));
        else
            iR = inv(class(k).R);
        end
        class(k).iR=iR;
        
        % Gradient
        gk = Upre*(c(:,k)-y(:,k))-class(k).R0*(class(k).w-class(k).w0);
        
        wold = class(k).w;
        dw = iR *gk;
        
        Jnew=-Inf;
        alpha=alpha_init;
        for i=1:10,
            class(k).w = wold + alpha*dw;
            bsr.class=class;
            if weighting
                Jnew = bsr_logjoint(bsr,U,c,gamma);
            else
                Jnew = bsr_logjoint(bsr,U,c);
            end
            if Jnew < Jlast
                alpha=alpha/2;
%                 if verbose
%                     fprintf('Halving step size ...\n');
%                 end
            else
                break
            end
        end
        if Jnew < Jlast,
            class(k).w=wold;
            break;
        end
        Jlast=Jnew;
    end
end

bsr.class=class;
bsr.labels=labels;



