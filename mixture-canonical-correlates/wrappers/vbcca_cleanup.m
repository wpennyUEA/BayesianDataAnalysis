function [cca_clean,latent_vars] = vbcca_cleanup (cca,prt)
% Clean up CCA model
% FORMAT [cca_clean,latent_vars] = vbcca_cleanup (cca,prt)
%
% cca           original model
% prt           remove clusters with pi(k) < prt
%
% cca_clean     cleaned up model
% latent_vars   vector containing numbers of latent variables in 
%               remaining clusters
%
% Remove empty (or sparse) clusters
% Remove empty columns of factor loading matrices

if nargin < 2, prt = 0.01; end

K = length(cca.W1);

cca_clean.Fhist = cca.Fhist;
cca_clean.F = cca.F;
cca_clean.Fdecomp = cca.Fdecomp;

m=1;
for k=1:K,
    %if sum(sum(abs(cca.W1{k}))) > 0 & sum(sum(abs(cca.W2{k}))) > 0

    W2 = cca.W2{k};
    W1 = cca.W1{k};
    p = size(W2,2);

    % Remove any columns of W2 that have shrunk towards zero
    % i.e. latent variable effectively removed
    rem_col = [];
    for jj = 1:p,
        if norm(W2(:,jj)) < 1e-3
            rem_col = [rem_col, jj];
        end
    end
    if length(rem_col) > 0
        W2(:,rem_col) = [];
        W1(:,rem_col) = [];
    end
     
    if (length(rem_col) < p) & (cca.pi(k) > prt)
        latent_vars(m) = size(W1,2);
        cca_clean.W1{m} = W1;
        cca_clean.W2{m} = W2;
        cca_clean.mu1{m} = cca.mu1{k};
        cca_clean.mu2{m} = cca.mu2{k};
        cca_clean.psi1{m} = cca.psi1{k};
        cca_clean.psi2{m} = cca.psi2{k};
        cca_clean.C1{m} = cca.C1{k};
        cca_clean.C2{m} = cca.C2{k};
        cca_clean.pi(m) = cca.pi(k);
        m = m+1;
    end
end
cca_clean.pi = cca_clean.pi/sum(cca_clean.pi);
cca_clean.pi = cca_clean.pi(:);
cca_clean.M = length(cca_clean.W1);