function [C1hat,gamma,p1cond,org] = vbcca_conditional_1D (cca,S,c2,con)
% Plot 1D state density as a function of context c2/x2
% FORMAT [C1hat,gamma,p1cond,org] = vbcca_conditional_1D (cca,S,c2,con)
%
% cca       cca data structure from vbcca.m
% S         (original) space to plot over 
%           .range [x1min, x1max];
%           .msd [mu1 sd1];
% c2        c2 vector to condition on; c2=x2 if con not specified
% % con       (optional)
%           .Gamma1  [1 x d1] contrast matrix s.t. c1 = Gamma1*x1
%           .Gamma2  [n2 x d2] contrast matrix s.t. c2 = Gamma2*x2
%
% Outputs:
%
% gamma     cluster assignments
% C1hat     predictions of C1 (i.e. mean of predictive density)
% p1cond    density over c1 in org space
% org       original space

bins = 100;
norm_range = norm_space (S.range,S.msd);
S.xmin = norm_range(1,1);
S.xmax = norm_range(1,2);
S.dx = (S.xmax-S.xmin)/bins;
x = [S.xmin:S.dx:S.xmax];
org = org_space(x,S.msd);
Nx = length(x);

Ypred = x';
X1p = Ypred;

if nargin < 4
    x2 = c2;
    [p1cond,gamma] = vbcca_conditional (cca,X1p,x2);
else
    if size(con.Gamma1,1) == 1
        [C1hat,gamma,p1cond,T] = vbcca_cond_subspace (cca,c2,con,X1p);
    else
        disp('Error in vbcca_conditional_1D: x1 subspace expected to be 1D');
        return
    end
end



 