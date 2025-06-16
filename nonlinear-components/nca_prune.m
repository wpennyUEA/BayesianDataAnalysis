function [M] = nca_prune (X, y, M)
% Bayesian Neighborhood Component Analysis with Pruning
% FORMAT [M] = nca_prune (X, y, M)
%
% X         [N x K] matrix of data
% y         [N x 1] Class labels
%           N is number of data points, K number of features
%
% M         .p              Initial Latent Dimension (default, p=K)
%           .opt            'LineSearch' (default) or 'FixedStepSize'
%           .lambda         Prior precision of feature weights (default=100)
%           .maxits         (default = 128)
%           .prune_its_min  (default = 8). Set to > .maxits to **turn pruning off**.
%           .verbose   0/1  (default = 1)
%
% Outputs are:
%
% M         .p          Final latent dimension
%           .A          [M.p x K] Encoding matrix
%           .c          Cost by iteration number
%           .cr         Correct Rate
%           .rep        Repetition number
%                
% [1] Goldberger et al. Neighborhood Component Analysis, NIPS 2004
%
% Remove rows of A which do not significantly reduce model evidence
% Uses batch (not online) learning
%
% To turn the pruning off, and just make use of the global shrinkage prior,
% set prune_its_min > maxits
%
% Learning rate for 'FixedStepSize' optimisation
alpha=0.1;

[N,K] = size(X);

if isfield(M,'p') p = M.p; else M.p=K; end
if isfield(M,'lambda') lambda=M.lambda; else lambda=100; end
if isfield(M,'verbose') verbose = M.verbose; else verbose=0; end
if isfield(M,'opt') opt = M.opt; else opt='LineSearch'; end
if isfield(M,'cost') cost = M.cost; else cost='gp'; end
if isfield(M,'maxits') maxits=M.maxits; else maxits=128; end

% Pruning opportunity every "gap" iterations after "min" iteration number
if isfield(M,'prune_its_min') prune_its_min=M.prune_its_min; else prune_its_min=8; end
prune_its_gap=8;

M.gradient=1;
M.cost='g';

% Initialisation
A = 0.1*randn(M.p,K);
grad=zeros(M.p,K);
M.lambda=lambda*ones(M.p,1);

% Max step size for linear search method
alpha_max =1;
nr_steps=0;
            
for it=1:maxits,
      
    M.gradient=1;
    [g,grad] = nca_org_cost_grad (X, y, M, A);
      
    % Add on prior contribution to gradient
    gp=g;
    for k=1:M.p,
        grad(k,:) = grad(k,:) - M.lambda(k)*A(k,:);
    end
    ng=norm(grad);
    dLdA = grad/ng;
    %dLdA = grad;
    
    oldA = A;
    switch opt,
        case 'FixedStepSize',
            A = A + alpha*dLdA;
            
        case 'LineSearch',
            [alpha, E] = fminbnd(@(alpha) tune_alpha(alpha,X,y,M,A,dLdA),0,alpha_max);
            A = A + alpha*dLdA;
            
        otherwise
            disp('Error in nca_prune.m: Unknown optimisation method');
            return
    end
    
    M.gradient=0;
    g = nca_org_cost_grad (X, y, M, A);
    lp = nca_log_prior (M.lambda,A);
    c(it) = g + sum(lp);
    
    if isnan(c(it))
        keyboard
    end
    
    if it > 1 & (c(it) < c(it-1))
        A = oldA;
        c(it) = c(it-1);
        normal_step=0;
        nr_steps=nr_steps+1;
        alpha_max = alpha_max/2;
        if verbose
            disp(sprintf('It %d Reduce max step size, Log Joint =%1.2f',it,c(it)));
        end
    elseif verbose
        normal_step=1;
        nr_steps=0;
        disp(sprintf('It %d Step Size = %1.4f, Log Joint =%1.2f',it,alpha,c(it)));
    end
    
    if it > 8
        small_prop_change = abs((c(it)-c(it-1))/c(it)) < 0.001;
        max_nr_steps = nr_steps > 3;
        if small_prop_change | max_nr_steps
            break;
        end
    end
    
    if (it > prune_its_min) 
        mm=0;
        if (mm==0)
            % Pruning opportunity
            M.verbose=0;
            removed=1;
            while (removed==1) & (M.p > 1)
                [A,removed] = nca_remove_row (X,y,M,A);
                M.p = size(A,1);
            end
        end
    end
end

if (it > prune_its_min) & (M.p > 1)
    % More pruning at end of optimisation
    M.verbose=0;
    removed=1;
    while removed==1
        [A,removed] = nca_remove_row (X,y,M,A);
        M.p = size(A,1);
    end
end

Mtmp=M;
Mtmp.cost='f';
total_correct = nca_org_cost_grad (X, y, Mtmp, A);
correct_rate=total_correct/N;

if verbose
    disp(' ');
    disp(sprintf('Correct Rate = %1.2f',correct_rate));
    disp(' ');
end

M.A = A;
M.p = size(A,1);
M.cr = correct_rate;


