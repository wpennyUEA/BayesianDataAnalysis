function [J,y] = blr_logjoint (w,M,X,t)
% Log Joint for Bayesian logistic regression
% FORMAT [J,y] = blr_logjoint (w,M,X,t)
%
% w     parameter vector
% M     priors
% X     input features
% t     binary class labels
%
% J     log joint
% y     sigmoid outputs

a = w'*X'; % activations
y = 1./(1+exp(-a)); % sigmoid output
y = y(:);

L = sum(t.*log(y)+(1-t).*log(1-y));

e=w-M.m0;
J = L -0.5*e'*M.R0*e;

