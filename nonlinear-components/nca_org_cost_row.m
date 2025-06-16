function [c] = nca_org_cost_row (a,X,y,M,A,k)
% Return cost as a function of kth row of A
% FORMAT [c] = nca_org_cost_row (a,X,y,M,A,k)
%
% a             what A(k,:) will be
% X,y,m,A       See nca_ord_cost_grad
% k             kth row 
%
% c             cost

M.gradient=0;
M.cost='g';
A(k,:)=a;
c = nca_org_cost_grad (X, y, M, A);