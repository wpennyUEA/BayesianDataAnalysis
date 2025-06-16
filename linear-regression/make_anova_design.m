function [X,y,group] = make_anova_design (group,Xc)
% Make design matrix and data vector for linear model
% FORMAT [X,y,group] = make_anova_design (group,Xc)
%
% group(g).x        variables
% group(g).name     name
% Xc                [N x C] matrix containing confounds to regress-out/control-for
%                   where C is number of confounding variables
%                   N is total number of data points
%
% X                 design matrix
% y                 data vector

G=length(group);

y=[];gr=[];
X=[];
Nchk=0;
for g=1:G,
    x=group(g).x(:);
    group(g).m=mean(x);
    N=length(x);
    Nchk=Nchk+N;
    group(g).sem=std(x)/sqrt(N-1);
    y=[y;x];
    gr=[gr;g*ones(N,1)];
    X=blkdiag(X,ones(N,1)); 
end
if nargin > 1
    % Add confounds if specified
    [N,C] = size(Xc);
    if ~(N==Nchk)
        disp('Error in make_anova_design.m');
        disp(sprintf('Rows in Xc (=%d) not equal to total number of datapoints (=%d)',N,Nchk));
        keyboard
        return
    end
    X = [X, Xc];
end