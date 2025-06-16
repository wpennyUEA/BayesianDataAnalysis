function [C1hat,gamma] = vbcca_conditional_2D_new (cca,S,c2,con)
% Plot 2D state density as a function of context c2/x2
% FORMAT [C1hat,gamma] = vbcca_conditional_2D_new (cca,S,c2,con)
%
% cca       cca data structure from vbcca.m
% S         (original) space to plot over 
%           .range [x1min x1max; x2min x2max];
%           .msd [mu1 sd1 ; mu2 sd2];
% c2        c2 vector to condition on; c2=x2 if con not specified
% con       (optional)
%           .Gamma1  [2 x d1] contrast matrix s.t. c1 = Gamma1*x1
%           .Gamma2  [n2 x d2] contrast matrix s.t. c2 = Gamma2*x2
%
% Outputs:
%
% gamma     cluster assignments
% C1hat     predictions of C1 (i.e. mean of predictive density)

bins = 30;
norm_range = norm_space (S.range,S.msd);
S.xmin = norm_range(2,1);
S.xmax = norm_range(2,2);
S.ymin = norm_range(1,1);
S.ymax = norm_range(1,2);
S.dx = (S.xmax-S.xmin)/bins;
S.dy = (S.ymax-S.ymin)/bins;

x = [S.xmin:S.dx:S.xmax];
y = [S.ymin:S.dy:S.ymax];

org = org_space([x;y],S.msd);
    
Nx = length(x);
Ny = length(y);

k = 1;
for i=1:Nx,
    for j=1:Ny,
        Ypred(:,k)= [x(i) y(j)]';
        k = k+1;
    end
end

X1p = Ypred;

if nargin < 4
    x2 = c2;
    [p1cond,gamma] = vbcca_conditional (cca,X1p,x2);
else
    if size(con.Gamma1,1) == 2
        [C1hat,gamma,p1cond,T] = vbcca_cond_subspace (cca,c2,con,X1p);
    else
        disp('Error in vbcca_conditional_2D: x1 subspace expected to be 2D');
        return
    end
end

k = 1;
for i=1:Nx,
    for j=1:Ny,
        D1(j,i)=p1cond(k);
        k = k+1;
    end
end

imagesc(org(1,:),org(2,:),D1);
axis xy
%title(sprintf('p(x1|x2=[%1.1f, %1.1f])',x2(1),x2(2)));
% disp(sprintf('Soft cluster assignments:'));
% disp(gamma);


 