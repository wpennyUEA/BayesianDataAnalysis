function [p,m,v] = mygampdf (x,b,c)
% Format PDF of gamma distribution
% FORMAT [p,m,v] = mygampdf (x,b,c)
%
% x     variate
% b 
% c
%
% 

L = (c-1)*log(x)-x/b-c*log(b)-gammaln(c);
p = exp(L);
m = b*c;
v = b^2*c;

