function [p1,gamma] = vbcca_conditional (cca,X1p,x2)
% Compute conditional density, p(x1|x2)
% FORMAT [p1,gamma] = vbcca_conditional (cca,X1p,x2)
%
% cca       data structure returned by vbcca.m
% X1p       [d1 x N] matrix of X1 data points over which to compute density
% x2        [d1 x 1]  data vector to condition on
%
% p1        [1 x N] vector of probability densities for X1p
% gamma     [M x 1] vector of cluster responsibilities

p=[];

for k=1:cca.M,
    
    % p(z^k|x_2,m=k)
    W2=cca.W2{k};
    p = size(W2,2);
    iC = cca.psi2{k};
    T1 = W2'*iC*W2;
    mu = inv(T1)*W2'*iC*(x2-cca.mu2{k});
    C = inv(T1+eye(p));
    
    % p(x1|x2,m=k)
    W1 = cca.W1{k};
    r1 = W1*mu + cca.mu1{k};
    C1 = cca.C1{k};
    Sigma1 = W1*C*W1'+ C1;
    Sigma1 = 0.5*(Sigma1+Sigma1'); % ensure symmetric
    pdf1(k,:) = spm_mvNpdf (X1p,r1,Sigma1);
    
    % p(x2|m=k)
    mu2 = cca.mu2{k};
    C2 = cca.W2{k}*cca.W2{k}' + cca.C2{k};
    C2 = 0.5*(C2+C2'); % ensure symmetric
    pdf2(k,:) = spm_mvNpdf (x2,mu2,C2);
       
end

% p(X1p|x2)
u = pdf2.*cca.pi;
gamma = u/sum(u);  % gating, p(m=k|x2)
p1 = gamma'*pdf1;
