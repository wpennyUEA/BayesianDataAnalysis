function [cca] = vbcca (x1,x2,k,M,options)
% Mixtures of Canonical Correlate Analysers
% FORMAT [cca] = vbcca (x1,x2,k,M,options)
%
% Inputs:
% x1        [d1 x N] data matrix
% x2        [d2 x N] data matrix
% k         maximum latent dimension
% M         number of clusters (default, M=1)
% options   .maxIter  maximum number of iterations (default = 128)
%           .tol      tolerance for convergence (default = 10^(-4))
%           .null     set to 1 for null model, 0 otherwise (default)
%           .retries  number of KNN (initialisation) retries
%           .indices  if using vbcca_index_init.m
%
% Outputs:
% cca   .W1 Factor loading matrix for x1   
%       .W2 Factor loading matrix for x2
%       .mu1 mean for x1
%       .mu2 mean for x2
%       .psi1 obs noise precision for x1
%       .psi2 obs noise precision for x2
%       .C1 obs noise covariance for x1
%       .C2 obs noise covariance for x2
%       .pi [M x 1] vector of prior cluster probs
%       .F  approx to model evidence
%       .Fhist history thereof
%       .cca_model  internal data structure
%
% Set M=1 for Bayesian CCA (default)
%
% This function calls code by [1]
%
% [1] J Viinikanoja, A Klami and S Kaski, 
% Variational Bayesian Mixture of Robust CCA models, ECML PKDD, 2010
%
% If options.null=1 the factor loading matrices are fixed at (close to) 0
% by setting a high prior precision (about 0). This permits computation
% of Bayes Factors relative to a null model.

try maxIter = options.maxIter; catch maxIter = 128; end
try tol = options.tol; catch tol = 10^(-4); end
try null_model = options.null; catch null_model = 0; end
if nargin < 4 | isempty(M), M=1; end
if M==1, 
    retries=1; 
else 
    try retries = options.retries; catch retries = 25; end
end

[d1,N] = size(x1);
[d2,N] = size(x2);

%kmax = min([d1,d2])-1;
kmax = min([d1,d2]);
if nargin < 3, k=kmax; end
if k > kmax,
    disp(sprintf('Latent dimension set to k=%d',k));
    k = kmax;
end

% Initialise model with normal latents (NOT t-distributed)
nModel = vbmcca(M,d1,d2,k);
if null_model
    mean_prior_precision = 1000;
    var_prior_precision = 1;
    nModel.ard_b = mean_prior_precision/var_prior_precision;
    nModel.ard_a = nModel.ard_b*mean_prior_precision;
end

nModel.normalLatents = 1;

if isfield(options,'indices')
    nModel = vbcca_index_init (nModel, x1, x2, options.indices);
else
    nModel = nModel.initWithKMeans(x1,x2,150,retries);
end


% Learning
%---------------------------------------------
learningParameters = vbmcca.initLearningParameters(maxIter, tol);
[cca_model, nEnergies, Fdecomp] = nModel.learn(x1, x2, learningParameters); 

cca.cca_model = cca_model;
cca.Fhist = nEnergies;
cca.F = nEnergies(end);
cca.Fdecomp = Fdecomp;
cca.M = M;
cca.pi = cca.cca_model.pi;

% Check this is the same as cca.F
% [1 1 -ones(1,11)]*Fdecomp.term(:)

for m = 1:M,
    cca.W1{m} = cca_model.eX_W(:,:,m);
    cca.W2{m} = cca_model.eY_W(:,:,m);
    cca.mu1{m} = cca_model.eX_mu(:,m);
    cca.mu2{m} = cca_model.eY_mu(:,m);
    cca.psi1{m} = cca_model.eX_psi(:,:,m); % obs noise precision 
    cca.psi2{m} = cca_model.eY_psi(:,:,m); % obs noise precision 
    cca.C1{m} = inv(cca.psi1{m}); 
    cca.C2{m} = inv(cca.psi2{m}); 
end



