function [H,M] = get_corr_entropy (Y,ind)
% Get correlation entropy using selected data only
% FORMAT [H,M] = get_corr_entropy (Y,ind)
%
% Y     [D x T] data matrix
% ind   [S x 1] vector of selected indices/time points
%       default is all indices
%
% H     Correlation entropy
% M     Optimal number of PCA components

if nargin < 2 | isempty(ind);
    Y=Y;
else
    Y=Y(:,ind);
end

C = cov(Y');
riC = diag(sqrt(1./diag(C)));
R = riC*C*riC;
H = gaussian_entropy (R);

if nargout > 1
    N = size(Y,2);
    M = spm_pca_order(C,N);
end