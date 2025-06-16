function [mlm] = spm_mlm_bayes (y,x,options)
% Bayesian Multivariate Linear Modelling
% FORMAT [mlm] = spm_mlm_bayes (y,x,options)
%
% MLM: y = x W + e
%
% y           N-by-d data matrix (dependent variables)
% x           N-by-p design matrix (independent variables)
% options       
% .pr          Shrinkage prior on MLM coefficients:
%             'input' (default), 'output' or 'global'
%
%             For 'input', coeffs of each independent variable
%             ie rows of W, share same prior precision. This 
%             allows some inputs to be more relevant than others.
%
%             For 'output', cols of W share same prior precision.
%             This allows some outputs to be more relevant.
%
%             For 'global' there is a single prior precision which
%             corresponds to a Bayesian equivalent of ridge regression.
%
%             For 'spatial-manova' there is a single regularisation parameter
%             that controls the magnitude/smoothness of condition means.
%             Observations (j=1..d) are assumed to be spatially related as 
%             defined by (an unscaled) precision matrix, S. See
%             bayes_manova1.m for further details.
%             
%
% .verbose     1 to print out iteration details, 0 otherwise (default=0)
% .ml_only     set to 1 to only compute ML solution. Default is zero
% .alpha.mean   Mean prior precision of coeffs
% .alpha.std    STD prior precision of coeffs
% .pseudo      1 to use pinv instead of inv (default 0). 
%              Useful for large data matrices.
% .S           Unscaled spatial precision matrix
%              
%
% The returned data structure mlm contains the following fields
%
% .wmean      Bayes estimate of [p x d] regression coefficient matrix
% .wsd        [p x d] posterior standard deviations of reg coeffs
% .wml        Maximum Likelihood regression coefficient matrix
% .wcov       [pd x pd] posterior covariance of regression coeffs
% .lambda     [d x d] observation noise precision matrix
% .acc        Model Accuracy
% .acc_uni    Univariate contribution to model accuracy
% .acc_multi  Multivariate contribution to model accuracy
% .complexity Model complexity
% .fm         Negative (variational) free energy of model
%             (this is approx to Bayesian model evidence)
% .bic        Bayesian Information Criterion
% .iterations Number of iterations during optimisation
% .prior      Details of regression coeff prior
%             .group(j).mean_alpha:
%             Estimated prior precision of jth parameter group.
%             For 'input' prior this is jth row of W. 
%             For 'output' prior this is jth column of W.
%
% Adapted from spm_mar.m by allowing user-specified data matrix, X, instead
% of fixing X to be a time lagged version of Y

pr_default = 'input';
if nargin < 3, 
    options=[];
    verbose=0;
    ml_only=0;
    pr = pr_default;
else
    if isfield(options,'pr'), pr=options.pr; else pr = pr_default; end
    if isfield(options,'verbose'), verbose=options.verbose; else verbose = 0; end
    if isfield(options,'ml_only'), ml_only=options.ml_only; else ml_only = 0; end
    if isfield(options,'alpha'), alpha=options.alpha; else alpha = []; end
    if isfield(options,'pseudo'), pseudo=options.pseudo; else pseudo = 0; end
end

d=size(y,2);    % Dimension of each dependent variable
N=size(y,1);    % Number of data points
p=size(x,2);    % Dimension of independent variable

if d >= N | p >= N
    disp('Error in spm_mlm_bayes.m: too few data points');
    return
end

k=p*d;
switch lower(pr)
    case {'input','spatial-manova'}
        % Separate precision for each row of W ie for each input
        prior.type='input';
        prior.groups=p;
        % get indices of each 'parameter group' ie rows
        vec_ind=[1:k]';
        ind=reshape(vec_ind,p,d);
        g=zeros(1,k);
        for j=1:p,
            gj=g; gj(ind(j,:))=1;
            group(j).index=gj;
        end
        
    case 'output',
        % Separate precision for each column of W ie for each output
        prior.type='output';
        prior.groups=d;
        % get indices of each 'parameter group' ie cols
        vec_ind=[1:k]';
        ind=reshape(vec_ind,p,d);
        g=zeros(1,k);
        for j=1:d,
            gj=g; gj(ind(:,j))=1;
            group(j).index=gj;
        end
        
    otherwise
        % Global prior
        prior.type='global';
        prior.groups=1;
        group(1).index=ones(1,k);
end

% Get both pseudo-inverse and approx to inv(x'*x) efficiently
[ux,dx,vx]=svd(x);
kk=size(x,2);
if kk==1
    xp=pinv(x);
    inv_xtx=1/(x'*x);
else
    ddx=diag(dx);
    svd_tol=max(ddx)*eps*kk;
    rank_X=sum(ddx > svd_tol);
    ddxm=diag(ones(rank_X,1)./ddx(1:rank_X));
    ddxm2=diag(ones(rank_X,1)./(ddx(1:rank_X).^2));
    xp=vx(:,1:rank_X)*ddxm*ux(:,1:rank_X)';  % Pseudo-inverse
    inv_xtx=vx(:,1:rank_X)*ddxm2*vx(:,1:rank_X)'; % approx to inv(x'*x)
end

% Compute terms that will be used many times
xtx=x'*x;
yty=y'*y;
xty=x'*y;
vec_xty=xty(:);

% Get maximum likelihood solution
w_ml = xp*y;
y_pred = x*w_ml;
H = y_pred'*y_pred;
if pseudo
    [tmp1,a,tmp2]=svd(H*pinv(yty),"econ");
else
    [tmp1,a,tmp2]=svd(H*inv(yty),"econ");
end

%[v_tmp,a]=eig(H*inv(yty));
mlm.stats.a = diag(a);
mlm.stats.pillai = sum(mlm.stats.a);
mlm.stats.r2 = mlm.stats.pillai/d;
if ~isreal(mlm.stats.a)
    keyboard
end

if ml_only
    mlm.wml=w_ml;
    return
end

e=y-y_pred;
noise_cov=(e'*e)/N;
sigma_ml=kron(noise_cov,inv_xtx);

% Priors on alpha(s)
if isempty(alpha)
    b_alpha_prior=1000;
    c_alpha_prior=0.001;
    mean_alpha_prior=b_alpha_prior*c_alpha_prior;
else
    var = alpha.std^2;
    b_alpha_prior = var/alpha.mean;
    c_alpha_prior = alpha.mean/b_alpha_prior;
    mean_alpha_prior = alpha.mean;
end

% Create spatial precision matrix specific to condition j=1..p
if strcmp(lower(pr),'spatial-manova')
    for j=1:prior.groups,
        S{j} = zeros(k,k);
        % ind(j,:)
        for aa = 1:d
            for bb = 1:d
                i1 = ind(j,aa);
                i2 = ind(j,bb);
                S{j}(i1,i2) = options.S(aa,bb);
            end
        end
    end
end


% Initialise 
w_mean=w_ml;
w_cov=sigma_ml;

max_iters=32;
w=zeros(p,d);
%tol=0.0001;
tol=1e-6;
for it=1:max_iters,
    
    %------------------------------------------------------------------
    % Update weight precisions
    if strcmp(lower(pr),'spatial-manova')
        b_alpha = 0;
        for j=1:prior.groups,
            b_alpha = 0.5*w_mean(:)'*S{j}*w_mean(:) + 0.5*trace(w_cov*S{j});
        end
        b_alpha = b_alpha +(1/b_alpha_prior);
        group(1).b_alpha = 1/b_alpha;
        group(1).c_alpha=0.5*prior.groups + c_alpha_prior;
        group(1).mean_alpha=group(1).b_alpha*group(1).c_alpha;
    else
        for j=1:prior.groups,
            Ij=diag(group(j).index);
            kj=sum(group(j).index);
            E_wj=0.5*w_mean(:)'*Ij*w_mean(:);
            b_alpha=E_wj+0.5*trace(Ij*w_cov*Ij)+(1/b_alpha_prior);
            group(j).b_alpha=1/b_alpha;
            group(j).c_alpha=0.5*kj+c_alpha_prior;
            group(j).mean_alpha=group(j).b_alpha*group(j).c_alpha;
            group(j).E_w=E_wj;
        end
    end

    yy_pred=x*w_mean;
    ee=y-yy_pred;
    E_d_av=ee'*ee;
    
    Omega = get_omega (p,d,w_cov,xtx);
    E_d_av=E_d_av+Omega;
       
    %------------------------------------------------------------------
    % Update noise precision posterior
    B=E_d_av;
    a=N;
    if pseudo
        mean_lambda=a*pinv(B);
    else
        mean_lambda=a*inv(B);
    end
    
    ilambda=(1/a)*B;
    
    % prior_cov=zeros(k,k);
    % if strcmp(lower(pr),'spatial-manova')
    %     for j=1:prior.groups,
    %         prior_cov=prior_cov+(1/group(j).mean_alpha)*Ij;
    %     end
    %     prior_cov = inv(prior_prec);
    % else
    %     for j=1:prior.groups,
    %         Ij=diag(group(j).index);
    %         prior_cov=prior_cov+(1/group(j).mean_alpha)*Ij;
    %     end
    % end

    %------------------------------------------------------------------
    % Convergence criterion
    old_w=w;
    w=w_mean;
    
    if (it<=10)
        w=w_mean;
    else
        change=norm(w(:)-old_w(:))/k;
        if verbose
            disp(sprintf('Iteration %d Delta_w=%1.6f',it,change));
        end
        if change < tol
            break;
        end
    end
   
    %------------------------------------------------------------------
    % Update weight posterior
    data_precision=kron(mean_lambda,xtx);
    prior_prec=zeros(k,k);

    if strcmp(lower(pr),'spatial-manova')
        for j=1:prior.groups,
            prior_prec=prior_prec+group(1).mean_alpha*S{j};
        end
    else
        for j=1:prior.groups,
            Ij=diag(group(j).index);
            prior_prec=prior_prec+group(j).mean_alpha*Ij;
        end
    end
    prior_cov = inv(prior_prec);

    if pseudo
        w_cov=pinv(data_precision+prior_prec);
    else
        w_cov=inv(data_precision+prior_prec);
    end

    vec_w_mean=w_cov*data_precision*w_ml(:);
    
    w_mean=reshape(vec_w_mean,p,d)
    
end

% Compute Negative Free Energy
kl_alpha=0;
if strcmp(lower(pr),'spatial-manova')
    kl_alpha=spm_kl_gamma(group(1).b_alpha,group(1).c_alpha,b_alpha_prior,c_alpha_prior);
else
    for j=1:prior.groups,
        kl_alpha=kl_alpha+spm_kl_gamma(group(j).b_alpha,group(j).c_alpha,b_alpha_prior,c_alpha_prior);
    end
end

kl_weights=spm_kl_eig_normal(w_mean(:),w_cov,prior_cov);
lga=spm_lg_gamma(d,0.5*a);
acc=-0.5*N*log(det(B));
f_m=acc-kl_weights-kl_alpha+lga;
if verbose
    disp('Contributions to Negative Free Energy');
    disp(sprintf('Acc=%1.3f KL_w=%1.3f KL_alpha=%1.3f LogGena=%1.3f Fm=%1.3f',acc,kl_weights,kl_alpha,lga,f_m));
end

% Get error bars on regression coefficients
w_sd=zeros(p,d);
post_var=diag(w_cov);
for dd=1:d,
    start=(dd-1)*p+1;
    stop=start+p-1;
    w_sd(:,dd)=sqrt(post_var(start:stop));
end

% Create mlm data structure
mlm.wmean=w_mean;
mlm.wsd=w_sd;
mlm.wml=w_ml;
mlm.wcov=w_cov;
mlm.lambda=mean_lambda;
mlm.acc=acc;
mlm.fm=f_m;
mlm.bic=-0.5*N*log(det(B))-0.5*k*log(N);
mlm.iterations=it;

mlm.kl_weights=kl_weights;
mlm.kl_alpha=kl_alpha;
mlm.lga=lga;
mlm.B = B;
mlm.prior=prior;
mlm.prior.group=group;
mlm.prior_cov = prior_cov;

% Univariate and multivariate contributions to model accuracy
Bdiag=diag(diag(B));
mlm.acc_uni = -0.5*N*log(det(Bdiag));
mlm.acc_multi = acc - mlm.acc_uni;

mlm.complexity = mlm.acc-mlm.fm;

function [Omega] = get_omega (p,d,w_cov,xtx)
% Get contribution to prediction error variance from w_cov 
% FORMAT [Omega] = spm_mlm_get_omega (p,d,w_cov,xtx)
%
% p         Number of independent variables
% d         Number of dependent variables
% w_cov     Uncertainty in MLM coefficients
% xtx       X'X where X is design matrix 
%
% Omega     Expected error variance from w_cov


Omega=zeros(d,d);
% Submatrix size - ie number of model params per dependent variable
s=p;

% Get upper diagonal elements
for di=1:d,
    for dj=di:d,
        istart=1+(di-1)*s;
        istop=istart+s-1;
        jstart=1+(dj-1)*s;
        jstop=jstart+s-1;
        w_cov_i_j=w_cov(istart:istop,jstart:jstop);
        Omega(di,dj)=trace(w_cov_i_j*xtx);
    end
end

% Get lower diagonal elements
Omega = Omega+Omega'-diag(diag(Omega));