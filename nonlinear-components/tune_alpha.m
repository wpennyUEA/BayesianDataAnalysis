function [E] = tune_alpha (alpha,X,y,M,A,dLdA)
% Evaluate negative log likelihood (error) 
% FORMAT [E] = tune_alpha (alpha,X,y,M,A,dLdA)
%
% X,y,M,A   See nca_prune.m
% alpha     Step size
% dLdA      Gradient
% 
% E         Negative log likelihood

M.gradient=0;
A = A + alpha*dLdA;
g = nca_org_cost_grad (X, y, M, A);
E = -g;