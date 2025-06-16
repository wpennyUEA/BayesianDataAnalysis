function [C2hat,gamma,p1,T] = vbcca_cond_context (cca,c1,con,C2p)
% Compute prediction and (optionally) conditional density in subspace, p(c2|c1)
% FORMAT [C2hat,gamma,p1,T] = vbcca_cond_context (cca,c1,con,C2p)
%
% cca       data structure returned by vbcca.m
% c1        [n1 x N1] data vectors to condition on
% con       .Gamma1  [n1 x d1] contrast matrix s.t. c1 = Gamma1*x1
%           .Gamma2  [n2 x d2] contrast matrix s.t. c2 = Gamma2*x2
%           where CCA was optimised for x1, x2 data space
% C2p       [n2 x N2] matrix of c2 data points over which to compute density
%           (optional)
%
% C2hat     [n2 x N1] vector of predictions
% gamma     [cca.M x N1] vector of cluster responsibilities
% p1        [1 x N2] vector of probability densities for X2p 
%           (only returned if C2p supplied)
% T         T{m} transformation matrix from c1 to c2 for mth cluster
%
% This function swaps around the x1,x2 parts of the cca model
% and then calls vbcca_cond_subspace

cca_swap = cca;
cca_swap.W2 = cca.W1;
cca_swap.W1 = cca.W2;
cca_swap.mu1 = cca.mu2;
cca_swap.mu2 = cca.mu1;
cca_swap.psi1 = cca.psi2;
cca_swap.psi2 = cca.psi1;
cca_swap.C1 = cca.C2;
cca_swap.C2 = cca.C1;
cca_swap.cca_model = [];
cca = cca_swap;

con_swap.Gamma1 = con.Gamma2;
con_swap.Gamma2 = con.Gamma1;
con = con_swap;

if nargin < 4
    [C2hat,gamma] = vbcca_cond_subspace (cca,c1,con);
else
    [C2hat,gamma,p1,T] = vbcca_cond_subspace (cca,c1,con,C2p);
end