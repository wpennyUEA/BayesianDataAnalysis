function [c,dcdA,pcorr,pr] = nca_org_cost_grad (X, y, M, A)
% Return cost and gradient for NCA over batch of samples
% FORMAT [c,dcdA,pcorr,pr] = nca_org_cost_grad (X, y, M, A)
%
% X         [N x K] matrix of inputs
% y         [N x 1] Outputs
% M         .p         Latent Dimension (default, p=K)
%           .lambda    Prior precision
%           .cost      'gp' (default), 'g', or 'f'
%           .gradient  1 to compute, 0 to not (default)
%
% c         cost
% dcdA      gradient
% pcorr     mean prob correct
% pr        neighbour probabilities
%
% [1] Goldberger et al. Neighborhood Component Analysis, NIPS 2004
%

dcdA=[];

[N,K] = size(X);
X=X';

if isfield(M,'p') p = M.p; else p=size(A,1); end
if isfield(M,'verbose') verbose = M.verbose; else verbose=1; end
if isfield(M,'cost') cost = M.cost; else cost='gp'; end
if isfield(M,'gradient') gradient = M.gradient; else gradient=0; end

lambda=M.lambda;

gr=zeros(K,K);
gr2=zeros(K,K);

for i=1:N,
    e=A*(X(:,i)-X);
    if p>1
        d=sum(e.^2);
    else
        d=e.^2;
    end
    %tmp=exp(-d);
    
    % Equation (1) in [1]
    %tmp(i)=0;
    %pr(i,:)=tmp/sum(tmp);
    
    % Subtract min distance from other data points
    % so that each point will always have at least one neighbour
    % (the closest one)
    stt=[1:N];
    stt(i)=[];
    dstt=d(stt);
    dstt=dstt-min(dstt);
    tmp_stt=exp(-dstt);
    
    tmp(stt)=tmp_stt;
    tmp(i)=0;
    pr(i,:)=tmp/sum(tmp);
   
    if sum(tmp)==0
        disp(sprintf('Error in nca_org_cost_grad.m: data point %d has no neighbours',i));
        keyboard
    end
    
    % Equation (2) in [1]
    ind=find(y==y(i));
    tmp2=pr(i,ind);
    prob(i)=sum(tmp2);
    
    if gradient
        % Compute inner terms of Eq (5) in [1]
        % ready for gradient computation
        term1=zeros(K,K);
        term2=zeros(K,K);
        for k=1:N,
            Xik=X(:,i)-X(:,k);
            dterm=pr(i,k)*Xik*Xik';
            term1=term1+dterm;
            if y(i)==y(k)
                term2=term2+dterm;
            end
        end
        
        % Inner part of Eq(5) in [1]
        gr=gr+prob(i)*term1-term2;
        
        % Inner part of Eq(7) in [1]
        if prob(i)>0
            gr2=gr2+term1-term2/prob(i);
        end
    end
    
end

% Expected number of points correctly classified
% Equation (3) in [1]
fcorr=sum(prob);
pcorr=fcorr/N;

% Equation (6) in [1]
g=sum(log(prob+eps));

gp=g;
for k=1:p,
    gp = gp - 0.5* lambda * A(k,:)*A(k,:)';
end

switch cost
    case 'f',
        c=fcorr;
    case 'g',
        c=g;
    case 'gp',
        c=gp;
end

if ~gradient
    return
else
    % Equation (7) in [1]
    dgdA = 2*A*gr2;
    switch cost
        case 'f',
            % Equation (5) in [1]
            dfdA = 2*A*gr;
            dcdA=dfdA;
        case 'g',
            dcdA=dgdA;
        case 'gp',
            % Prior
            dpdA = zeros(p,K);
            for k=1:p,
                dpdA(k,:)=dpdA(k,:)-lambda*A(k,:);
            end
            dcdA=dgdA+dpdA;
    end
end



