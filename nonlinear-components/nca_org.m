function [M] = nca_org (X, y, M)
% Neighbourhood Component Analysis - Maximum Likelihood optimisation
% FORMAT [M] = nca_org (X, y, M)
%
% X         [N x K] matrix of inputs
% y         [N x 1] Outputs
% M         .p         Latent Dimension (default, p=K)
%           .opt       'LineSearch' (default) or 'FixedStepSize'
%           .cost      'gp' (default) or 'f'
%           .racc      Retrain if accuracy < racc (default=0.65)
%           .verbose   0/1 (default = 1)
%
% M         .A      [M.p x K] Encoding matrix
%           .c      Cost by iteration number
%           .cr     Correct Rate
%           .rep    Repetition number
%                
% [1] Goldberger et al. Neighborhood Component Analysis, NIPS 2004
%
% In this full batch version we iterate through the whole data set
% multiple times as necessary until convergence.
%
% This algorithm can be susceptible to local maxima.
% Estimation restarts from different (random) initialisation
% if accuracy less then "racc". Final accuracy is "cr" achieved
% after "rep" repetitions (best value is returned).

% Learning rate for fixed step size option
alpha=0.1;

[N,K] = size(X);

if nargin < 3, M=[]; end
if isfield(M,'p') p = M.p; else p=K; end
if isfield(M,'verbose') verbose = M.verbose; else verbose=0; end
if isfield(M,'racc') racc = M.racc; else racc=0.65; end
if isfield(M,'cost') cost = M.cost; else cost='gp'; end
if isfield(M,'opt') opt = M.opt; else opt='LineSearch'; end
    
maxreps=8;
maxits=64;
M.gradient=1;
cr_best=-1;

M.lambda = 0;

for rep=1:maxreps,
    
    % Initialisation
    Arnd = 0.1*randn(p,K);
    if rep==1, 
        Ainit(:,:,rep) = Arnd; 
    else
        % Make initialisation orthogonal to previous one
        ip=Ainit(:,:,rep-1)*Arnd';
        Ainit(:,:,rep) = Arnd-ip*Arnd;
    end
    A=Ainit(:,:,rep);
    
    for it=1:maxits,
        
        [c(it),grad] = nca_org_cost_grad (X, y, M, A);
        ng=norm(grad);
        dLdA = grad/ng;
        
        switch opt,
            case 'FixedStepSize',
                A = A + alpha*dLdA;
                
            case 'LineSearch',
                alpha_max =1;
                [alpha, E] = fminbnd(@(alpha) tune_alpha(alpha,X,y,M,A,dLdA),0,alpha_max);
                A = A + alpha*dLdA;
                
            otherwise
                disp('Error in nca_prune.m: Unknown optimisation method');
                return
        end
        
        if verbose
            disp(sprintf('It %d Cost=%1.2f',it,c(it)));
        end
        if it > 8
            if abs((c(it)-c(it-1))/c(it)) < 0.001
                break;
            end
        end
    end
    
    Mtmp=M;
    Mtmp.cost='f';
    total_correct = nca_org_cost_grad (X, y, Mtmp, A);
    correct_rate=total_correct/N;
    
    if verbose
        disp(' ');
        disp(sprintf('Rep %d, Correct Rate = %1.2f',rep,correct_rate));
        disp(' ');
    end
    if correct_rate >= racc,
        A_best=A;
        cr_best=correct_rate;
        c_best=c;
        grad_best=grad;
        break; 
    else
        % Store best solution so far
        if correct_rate > cr_best
            A_best=A;
            cr_best=correct_rate;
            c_best=c;
            grad_best=grad;
        end
    end
end

M.A=A_best;
M.c=c_best;
M.rep=rep;
M.cr=cr_best;


