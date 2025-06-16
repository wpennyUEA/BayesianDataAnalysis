function [model,hier,D,F] = vbmfx (model,model_fit,hier,opt)
% Fit models to data from a group of subjects using Bayesian mixed-effects 
% FORMAT [model,hier,D,F] = vbmfx (model,model_fit,hier,opt)
%
% INPUTS:
%
% model{s}.mu   prior mean for subject s
%         .S    prior precision for subject s
%         .Y    data for subject s
%
% model_fit     name of routine to fit individual subjects' model
%               e.g. 'linear_fit'. This routine must take
%               model{s}.mu, .S, and .Y as input and return
%               .w  posterior mean
%               .R  posterior precision
%               .F  model evidence
%
% hier          option to set group level priors on 
%               mean via .mu0, .S0 and
%               precision via .a0, .b0
%
% opt           .max_init_subj_its
%               .max_subj_its
%               .max_group_its
%               .max_outer_loop_its 
%
%               See code below for default values
%
% OUTPUTS:
%
% model{s}.w    posterior mean for subject s
%         .R    posterior precision for subject s
%
% 
% hier          See vbmfx_group.m
%
% D             .indiv_model  Initial fitted models (prior to MFX)
%               .w_indiv(:,s) Initial Parameter estimates for subj s
%               .w_hier(:,s)  Final parameter estimates for sub s
%
% F             group model evidence
%
% This code is based on 
% [1] J. Daunizeau (2019) Variational Bayesian Modelling of Mixed-Effects, ArXiv:1903.09003

try max_init_subj_its=opt.max_init_subj_its; catch max_init_subj_its=16; end
try max_subj_its=opt.max_subj_its; catch max_subj_its=16; end
try max_group_its=opt.max_group_its; catch max_group_its=16; end
try max_outer_loop_its=opt.max_outer_loop_its; catch max_outer_loop_its=1; end

try hier=hier; catch hier=[]; end

Nsub = length(model);

% --------------------------------------------------------
% Fit individual subject models with initial priors e.g. shrinkage priors
sopt.maxits = max_init_subj_its;
for s=1:Nsub,
    M.pE = model{s}.mu;
    M.pC = inv(model{s}.S);
    model{s} = feval(model_fit,model{s},sopt);
    w_indiv(:,s) = model{s}.w;
end
indiv_model = model;

% --------------------------------------------------------
% Mixed Effects
sopt.maxits = max_subj_its;
gopt.maxits = max_group_its;

hier = vbmfx_group (model,hier,gopt);

for it = 1:max_outer_loop_its,
    disp(sprintf('MFX outer loop %d ...',it));
    for s=1:Nsub,
        % Fit individual subject models with hierarchical priors
        model{s}.mu = hier.mu;
        % model{s}.S = hier.S; % OLD ERRONEOUS CODE !!!
        model{s}.S = hier.Rexp;
        model{s} = feval(model_fit,model{s},sopt);
        w_hier(:,s) = model{s}.w;
    end
    % Update hierarchical priors
    hier = vbmfx_group (model,hier,gopt);
end

D.indiv_model = indiv_model;
D.w_indiv = w_indiv;
D.w_hier = w_hier;

% Compute Approx to Model Evidence
F = 0;
for s=1:Nsub,
    F = F + model{s}.F;
end
F = F - hier.kl_group_mean - hier.kl_group_precision;

