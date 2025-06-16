function [logBF,glm] = bayes_glm_anova1 (group, method, Xc, verbose)
% Bayesian GLM implementation of One-Way ANOVA
% FORMAT [logBF,glm] = bayes_glm_anova1 (group, method, Xc, verbose)
%
% group(g).x        variables
% group(g).name     name
% method            'G'     Multivariate Gaussian priors and posteriors
%                   'T'     Multivariate T priors and posteriors
% Xc                [N x C] matrix containing confounds to regress-out/control-for
%                   where C is number of confounding variables
%                   N is total number of data points
% verbose           default=0
%
% logBF             Log Bayes Factor in favour of alternative hypothesis
% glm               Linear Model parameters
%
% Single variance component assumed

if nargin < 3 | isempty(Xc)
    confounds=0;
else
    confounds=1;
end

if nargin < 4 | isempty(verbose)
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

glm.X = X;
[N,K] = size(X);

% Subtract data mean
my = mean(y);
vy = std(y)^2;
y = y-my;

switch method
    case 'T',
        % c0 must be larger than 1 for proper density
        % 2*c0 is the degrees of freedom of T-prior
        % (c0=1 for Cauchy)
        glm.c0=1.1;
        glm.B0=eye(K);
        glm.w0=zeros(K,1);
        
        % scale observation noise level to data.
        % this then changes scale of prior variance over regression
        % coefficients through T-prior
        glm.b0 = 1/(glm.c0*vy);

    case 'G',
        glm.pE=zeros(K,1);
        glm.Ce=vy*eye(N); % IID Error matrix
        % Scale prior variance by data variability
        glm.pC=vy*eye(K);
        
    otherwise
        disp('Unknown method in bayes_glm_anova1.m');
        return
end

% Single variance component assumed
glm = bayes_linear_estimate (glm,y,method);

logBF = bayes_linear_test (glm,c,method);

if verbose
    disp(' ');
    disp(sprintf('LogBF in favour of alternative = %1.2f', logBF));
end


