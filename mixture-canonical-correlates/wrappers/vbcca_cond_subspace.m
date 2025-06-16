function [C1hat,gamma,p1,T,std_dev] = vbcca_cond_subspace (cca,c2,con,C1p)
% Compute prediction and (optionally) conditional density in subspace, p(c1|c2)
% FORMAT [C1hat,gamma,p1,T,std_dev] = vbcca_cond_subspace (cca,c2,con,C1p)
%
% cca       data structure returned by vbcca.m
% c2        [n2 x N2] data vectors to condition on
% con       .Gamma1  [n1 x d1] contrast matrix s.t. c1 = Gamma1*x1
%           .Gamma2  [n2 x d2] contrast matrix s.t. c2 = Gamma2*x2
%           where CCA was optimised for x1, x2 data space
% C1p       [n1 x N1] matrix of c1 data points over which to compute density
%           (optional)
%
% C1hat     [n1 x N2] vector of predictions
% gamma     [cca.M x N2] vector of cluster responsibilities
% p1        [1 x N1] vector of probability densities for X1p 
%           (only returned if C1p supplied)
% T         T{m} transformation matrix from c2 to c1 for mth cluster
% std_dev   [n1 x k] matrix containing std_dev of predictions for kth
%           cluster

n1 = size(con.Gamma1,1);
[n2,N2] = size(c2);

if nargin > 3, output_density=1; else output_density=0; end
p1=[];

if N2==1
    pred = zeros(n1,cca.M);
else
    pred = zeros(n1,N2,cca.M);
end

for k=1:cca.M,
    
    W1 = cca.W1{k};
    W2 = cca.W2{k};

    % p(z^k|c_2,m=k)
    W2 = con.Gamma2*W2;
    
    mu2tilde = con.Gamma2*cca.mu2{k};
    C2tilde = con.Gamma2*cca.C2{k}*con.Gamma2';
    
    iC = inv(C2tilde);
    T1 = W2'*iC*W2;
    iT1 = inv(T1);
    
    % p(x1|x2,m=k)
    W1 = con.Gamma1*W1;

    % Transformation matrix
    T{k} = W1*iT1*W2'*iC;

    % mu = iT1*W2'*iC*(c2-mu2tilde);
    % C = inv(T1+eye(p));
    
    mu1tilde = con.Gamma1*cca.mu1{k};
    C1tilde = con.Gamma1*cca.C1{k}*con.Gamma1';
    
    r1 = T{k}*(c2-mu2tilde) + mu1tilde;

    if N2==1
        pred(:,k) = r1;
    else
        pred(:,:,k) = r1;
    end
    
    if output_density
        Sigma1 = W1*C*W1'+ C1tilde;
        Sigma1 = 0.5*(Sigma1+Sigma1'); % ensure symmetric
        std_dev(:,k) = sqrt(diag(Sigma1));
        pdf1(k,:) = spm_mvNpdf (C1p,r1,Sigma1);
    end
    
    % p(x2|m=k)
    C2 = W2*W2' + C2tilde;
    C2 = 0.5*(C2+C2'); % ensure symmetric
    pdf2(k,:) = spm_mvNpdf (c2,mu2tilde,C2);
end

% c1hat
u = pdf2.*cca.pi;
%gamma = u./(ones(n2,1)*sum(u,1));  % gating, p(m=k|c2)
gamma = u./(ones(cca.M,1)*sum(u,1));  % gating, p(m=k|c2)
if N2==1,
    C1hat = pred*gamma;
    %C1hat(1,:) = sum(gamma.*squeeze(pred(1,:,:)));
else
    if cca.M == 1
        C1hat = pred;
    else
        for d=1:n1,
            C1hat(d,:) = sum(gamma'.*squeeze(pred(d,:,:)),2);
        end
    end
end

if output_density
    % p(C1p|c2)
    p1 = gamma'*pdf1;
end






    