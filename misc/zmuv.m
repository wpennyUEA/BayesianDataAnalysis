function [X,m,s] = zmuv (X)
% Make all columns of X have zero mean and unit variance
% FORMAT [X,m,s] = zmuv (X)
%
% X     [N x P] matrix of original data
%
% X     [N x P] matrix of normalised data 
% m     [P x 1] vector of original means
% s     [P x 1] vector of original SDs

P=size(X,2);
for p=1:P,
    x=X(:,p);
    m(p)=mean(x);
    s(p)=std(x);
    X(:,p)=(x-m(p))/s(p);
end