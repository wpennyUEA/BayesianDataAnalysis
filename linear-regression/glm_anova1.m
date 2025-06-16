function [p,stats] = glm_anova1 (group, Xc, verbose)
% GLM implementation of One-Way ANOVA
% FORMAT [p,stats] = glm_anova1 (group, Xc, verbose)
%
% group(g).x        variables
% group(g).name     name
% Xc                [N x C] matrix containing confounds to regress-out/control-for
%                   where C is number of confounding variables
%                   N is total number of data points
% verbose           (default =0)
%
% p                 p-value
% stats             see glm_test_hypothesis.m plus
%                   .X      design matrix
%                   .y      data vector
%
% Single variance component assumed

if nargin < 2 | isempty(Xc)
    confounds=0;
else
    confounds=1;
end

if nargin < 3 | isempty(verbose)
    verbose=0;
end

if confounds
    C = size(Xc,2);
    [X,y,group] = make_anova_design (group,Xc);
else
    [X,y,group] = make_anova_design (group);
end

G=length(group);
if verbose
    disp(' ');
    for g=1:G,
        disp(sprintf('Mean %s = %1.3f, SEM = %1.3f',group(g).name,group(g).m,group(g).sem));
    end
end
    
Con = spm_make_contrasts(G);
c = Con(2).c; % Main effect of group is the 2nd contrast

if confounds
    % Pad contrast matrix with zeros for confounds
    [d,k]=size(c);
    c=[c,zeros(d,C)];
end

% Single variance component assumed.
% To change this we need to specify V and add as
% 4th argument to function below
stats = glm_test_hypothesis (X,y,c');

p = stats.p;
stats.X = X;
stats.y = y;

df1=round(stats.df(1));
df2=round(stats.df(2));

if verbose
    disp(' ');
    disp(sprintf('Main Effect: F(%d,%d)=%1.2f, p=%g', df1,df2,stats.F, p));
    disp(' ');
end

