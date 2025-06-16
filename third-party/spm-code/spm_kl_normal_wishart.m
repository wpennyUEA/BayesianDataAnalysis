function [dkl] = spm_kl_normal_wishart (Q,P)
% KL divergence between Normal-Wisharts
% FORMAT [dkl] = spm_kl_normal_wishart (Q,P)
%
% Q     .mu
%       .lambda
%       .a
%       .B
%
% P     same structure
%
% dkl   KL-Divergence between Q and P where
%       q(x,Gamma) = N(x;mu,lambda*Gamma) Wi(Gamma;a,B)
%       p(x,Gamma) has same form and
%       x is e.g. mean and lambda*Gamma is e.g. precision

% Expected KL over Normals
d=length(Q.mu);
bar_Gamma = Q.a*inv(Q.B);
e = Q.mu-P.mu;
e = e(:);
r = P.lambda/Q.lambda;
exp_KL_N = P.lambda*e'*bar_Gamma*e+0.5*d*(r-log(r)-1);

% KL over Wisharts
klW = spm_kl_wishart (Q.a,Q.B,P.a,P.B);

dkl = exp_KL_N + klW;

