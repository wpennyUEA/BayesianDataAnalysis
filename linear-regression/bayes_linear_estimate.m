function [glm] = bayes_linear_estimate (glm,y,method)
% Bayesian parameter estimation for Linear Model
% FORMAT [glm] = bayes_linear_estimate (glm,y,method)
%  
% glm       Input fields are
%
%           .X      [N x p] design matrix of independent variables
%
% y         [N x 1] vector containing dependent variable
%
% method    'G'     Multivariate Gaussian priors and posteriors
%           'T'     Multivariate T priors and posteriors
%
% glm       Output fields are
%
%           .Ep     Posterior Mean
%           .F      Log Model Evidence
%           .R2     Proportion of Variance Explained by Model
%
% See bayes_tlinear_estimate.m and bayes_glinear_estimate.m for more
% info on input and output fields

switch method
    case 'G',
        glm = bayes_glinear_estimate (glm,y);
        
    case 'T',
        glm = bayes_tlinear_estimate (glm,y);
        glm.Ep = glm.wN;
        
    otherwise
        disp('Unknown method in bayes_linear_estimate.m');
        return
end

var_y = std(y)^2;
e = y - glm.X*glm.Ep;
var_e = std(e)^2;
glm.R2 = (var_y-var_e)/var_y;

