function [logbf] = bayes_linear_test (glm,c,method,a)
% Bayesian hypothesis testing for Linear Model 
% FORMAT [logbf] = bayes_linear_test (glm,c,method,a)
%
% glm       output of bayes_linear_estimate.m
% c         [K x P] contrast matrix where P is number of regression
%           coefficients and K is dimension of contrast
% method    'G' or 'T'
% a         value of contrast under null hypothesis (zeros(K,1) by default)
%
% logbf     Log Bayes Factor in favour of alternative hypothesis
%
% Bayes Factors are computed using the Savage-Dickey approximation
%

[K,P] = size(c);
if nargin < 4 | isempty(a)
    a=zeros(K,1);
else
    Kchk=length(a);
    if ~(K==Kchk), disp('Contrast matrix and hypothesized value have inconsistent dimensions'); return; end
end

switch method
    case 'G',
        logbf = bayes_glinear_test (glm,c,a);
        
    case 'T',
        logbf = bayes_tlinear_test (glm,c,a);
        
    otherwise
        disp('Unknown method in bayes_linear_test.m');
        return
end