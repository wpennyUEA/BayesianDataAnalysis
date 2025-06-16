function [A,removed] = nca_remove_row (X,y,M,A) 
% Remove a row of NCA matrix using greedy search
% FORMAT [A,removed] = nca_remove_row (X,y,M,A) 

[p,K] = size(A);
if p==0 | K==0
    disp('Error in nca_remove_row: A has zero dimension');
    keyboard
end

M.gradient=0;
g = nca_org_cost_grad (X, y, M, A);
lp = nca_log_prior (M.lambda,A);
c = g + sum(lp);
    
for k=1:p,
    a = A(k,:);
    
    % Log Likelihood
    c0(k) = nca_org_cost_row (zeros(1,K),X,y,M,A,k);
    
    % Add on Log Prior
    lpk=lp;lpk(k)=[];
    c0(k)=c0(k)+sum(lpk);
    
    % Posterior precision matrix for kth row
    H = spm_diff ('nca_org_cost_row',a,X,y,M,A,k,[1 1]);
    iC(:,:,k) =-full(H)+M.lambda(k)*eye(K);
    
    ldC(k)=spm_logdet(squeeze(iC(:,:,k)));
end

% Log evidence of Full Model under Factorised Laplace
Fa = c-0.5*sum(ldC)+0.5*p*K*log(2*pi);

% Log evidence of row-reduced models
for k=1:p,
    Fp(k) = c0(k)-0.5*sum(ldC)+0.5*ldC(k)+0.5*(p-1)*K*log(2*pi);
end

% Log Bayes Factors
logbf=Fp-Fa;
[maxlogbf,ind]=max(logbf);
if maxlogbf > 3
    disp(sprintf('Removing row, LogBF=%1.2f',maxlogbf));
    A(ind,:)=[];
    M.p = M.p-1;
    removed=1;
else
    removed=0;
end

if M.verbose
    disp(' ');
    disp('Statistical tests for removing rows - reduced versus full:');
    disp(' ');
    disp('LogBF:');
    disp(logbf);
    disp('LogLF:');
    disp(c0-c);
end