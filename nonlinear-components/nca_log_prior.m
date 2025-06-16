function [lp] = nca_log_prior (lambda,A)
% Compute log prior for each row of A 
% FORMAT [lp] = nca_log_prior (lambda,A)

P=size(A,1);
for k=1:P,
    lp(k) = -0.5*lambda(k)*A(k,:)*A(k,:)';
end
