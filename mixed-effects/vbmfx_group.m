function [hier] = vbmfx_group (model,hier,opt)
% Estimate group level parameters
% FORMAT [hier] = vbmfx_group (model,hier,opt)
%
% INPUTS:
%
% model{k}.w    posterior mean for kth model
%         .R    posterior precision for kth model
%
% hier      .mu0, .S0
%           .a0, .b0 
%
% opt       .maxits    max number of iterations
%
% OUTPUTS:
%
% hier      .mu, .S
%           .a, .b
%
%           .Rexp       expected precision matrix
%
% p(w) = N(mu0,S0), q(w)=N(mu,S)
% p(r_p) = Ga(a0,b0), q(r_p)=Ga(ap,bp)
%
% New empirical prior for within-subject model fitting should have 
% mean .mu and precision .Rexp
%
% This code is based on 
% [1] J. Daunizeau (2019) Variational Bayesian Modelling of Mixed-Effects, ArXiv:1903.09003

K=length(model);
P=length(model{1}.w);

% Default priors on mean
try,
    hier.mu0=hier.mu0;
    hier.S0 =hier.S0;
catch
    hier.mu0=zeros(P,1);
    hier.S0=0.01*eye(P);
end

% Default priors on precision
try, 
    hier.b0 = hier.b0;
    hier.a0 = hier.a0;
catch
    mean_prec=10;
    var_prec=1000;
    b0 = mean_prec/var_prec;
    a0 = b0*mean_prec;
    hier.b0 = b0*ones(P,1);
    hier.a0 = a0*ones(P,1);
end

% Current Rexp value
try, 
    hier.Rexp = hier.Rexp;
catch
    hier.a = hier.a0;
    hier.b = hier.b0;
    hier.Rexp = diag(hier.a./hier.b);
end

sum_wk=zeros(P,1);
for k=1:K,
    sum_wk = sum_wk+model{k}.w;
    wmat(:,k) = model{k}.w;

    if any(diag(model{k}.R<0))
        disp('Error in vbmfx_group.m:');
        disp(sprintf('Model %d has negative posterior variance estimates',k));
        return
    end
    model{k}.iR = inv(model{k}.R);
end

for it=1:opt.maxits,
    
    % Update second-level mean
    hier.S = hier.S0+K*hier.Rexp;

    % For Debugging *************************************************
    if any(diag(hier.S)<0)
        keyboard;
    end
    hier.iS = inv(hier.S);
    
    if it>1, old_mu=hier.mu; end
    hier.mu = hier.iS*(hier.S0*hier.mu0+hier.Rexp*sum_wk);
    
    if it>1,
        dmu = hier.mu-old_mu;
        prop_delta = mean(abs(dmu))/mean(abs(old_mu));
        disp(sprintf('Iteration %d propchange = %1.4f',it,prop_delta));
        if prop_delta<0.001, break; end
    end
    
    % Update second-level precision
    for p=1:P,
        hier.a(p) = hier.a0(p)+0.5*K;
        bp = hier.b0(p);
        dbp = 0;
        for k=1:K,
            dbp = dbp+(model{k}.w(p)-hier.mu(p))^2+hier.iS(p,p)+model{k}.iR(p,p);
        end
        hier.b(p)=bp+0.5*dbp;
    end
    % For Debugging *************************************************
    if any(hier.b<0)
        keyboard;
    end
    hier.Rexp=diag(hier.a./hier.b);
end

% Compute KL-divergences
hier.kl_group_precision = 0;
for p=1:P,
    b = 1/hier.b(p); % spm_kl_gamma works with 1/b
    b0 = 1/hier.b0(p);
    c = hier.a(p);
    c0 = hier.a0(p);
    hier.kl_group_precision = hier.kl_group_precision + spm_kl_gamma (b,c,b0,c0);
end
C = inv(hier.S);
C0 = inv(hier.S0);
hier.kl_group_mean = spm_kl_normal (hier.mu,C,hier.mu0,C0);