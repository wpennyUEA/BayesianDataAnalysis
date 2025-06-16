function [model] = linear_fit (model,opt)
% Fit linear model to data
% FORMAT [model] = linear_fit (model,opt)
%
% INPUT:
%
% model.mu       prior mean 
%      .S        prior precision 
%      .Y        data
%
% opt           optimisation options (see vbmfx.m)
%
% OUTPUT:
%
% model.w        posterior mean 
% model.R        posterior precision 
% model.F        approximation to model evidence

M = model.M;
U = model.U;

M.pE = model.mu;
M.pC = inv(model.S);

[Ep,Cp,F] = linear_post(M,U,model.Y);
model.w = Ep;
model.R = inv(Cp);
model.F = F;