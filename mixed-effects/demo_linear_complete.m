
clear all
close all

% This script simulates hierarchical linear mixed effects data, and
% performs Bayesian parameter estimation at both a group and indivudla
% level using Variational Bayes

% Generate Data
sim.v = [2,1,1];    % True population level effects
sim.Lambda = diag([100 10 100]); sim.C = inv(sim.Lambda);   % Precision matrix for between-subject variability

% Within-subject parameters
sim.N = 5;  % Observations per subject
sd = 0.1;   % Noise SD for within-subject data
%sd = 1;
sim.lambda = 1/(sd^2);  % Precision of noise

% Number of parameters
P = length(sim.v);
%Prior mean and covariance for group-level effects (hierarchical prior)
%hier.mu0 = 2*ones(P,1);
hier.mu0 = zeros(P,1);  % Prior mean = 0
hier.S0 = 0.01*eye(P);  % Weak prior covariance

% Generate design matrix for linear model
[M,U] = linear_model (sim.N,sim.lambda);

% Number of subjects
Nsub=90;

% Simulate true individual parameters for each subject
w_tmp = spm_normrnd(sim.v,sim.C,Nsub);
Yall=[];    % Container for data across subjects

% Generate data for each subject based on their true parameters
for s=1:Nsub,
    model{s}.w_true = w_tmp(:,s);
    model{s}.Y = U.X*model{s}.w_true+sqrt(M.Ce)*randn(sim.N,1);
    Yall = [Yall,model{s}.Y];
end

% Plot individual subject data time courses
plot_data=0;
if plot_data
    rN = ceil(sqrt(Nsub));
    for s=1:Nsub,
        subplot(rN,rN,s);
        plot(U.t,model{s}.Y,'x');
        ylim([0 max(max(Yall))]);
        grid on
    end
end

% Prep model structures for Variational Bayes fitting
for s=1:Nsub,
    model{s}.M = M; % Model structure
    model{s}.U = U; % Input design
    model{s}.mu = M.pE; % Prior mean for parameters
    model{s}.S = inv(M.pC); % Prior covariance
end

% Use tight prior precision
tight_prior=0;
if tight_prior
    % Set tight group precision prior
    sigma0 = 1.68;
    mean_prec = ones(P,1)*(1/sigma0)^2;
    var_prec = 0.01*ones(P,1);
    hier.b0 = mean_prec./var_prec;
    hier.a0 = hier.b0.*mean_prec;
end

% Perform Variational Bayes mixed effects modelling
[model,hier,D,F] = vbmfx(model,'linear_fit',hier);

% Calculate average individual-level parameter estimates across subjects
indiv_mu = mean(D.w_indiv')';

% Plot group-level parameter estimates
h=figure;
set(h,'Name','Group Level');
P = length(sim.v);
rN = ceil(sqrt(P));
for p=1:P,
    theta(p,:) = [sim.v(p),hier.mu(p),indiv_mu(p)];
end
bar(theta);
legend('True','Hierarchical','Shrinkage');
grid on
ylim([0.75*min(sim.v) 1.25*max(sim.v)]);
xlabel('Parameter');

% Plot subject-level parameter estimates against true values
h=figure;
set(h,'Name','Subject Level');
for p=1:P,
    subplot(rN,rN,p);
    x = sort(w_tmp(p,:));
    plot(x,x,'k-');
    hold on
    plot(w_tmp(p,:),D.w_hier(p,:),'ro');    % Hierarchical estimates
    plot(w_tmp(p,:),D.w_indiv(p,:),'bo');   % Shrinkage estimates
    grid on
    legend('True','Hierarchical','Shrinkage');
end

% Calculate sum of squared errors for parameter estimates
disp(' ');
disp('Errors in within-subject parameter estimates:');
for p=1:P,
    e=D.w_hier(p,:)-w_tmp(p,:);
    Eh = sum(e.^2);
    e=D.w_indiv(p,:)-w_tmp(p,:);
    Es = sum(e.^2);
    disp(sprintf('Parameter %d SSE: Hier = %1.3f, Shrinkage = %1.3f',p,Eh,Es));
    sse(p,:) = [Eh,Es];
end

% Plot sum of squared errors for hierarchical vs shrinkage estimates
h = figure;
set(h,'Name','Within Subject Error');
bar(sse);
legend('Hierarchical','Shrinkage');
grid on
xlabel('Parameter');
ylabel('Sum Squared Error');

% Plot prior and posterior group precision estimates
disp(' ');
disp('Prior and Posterior Group Precisions:');
vbmfx_plot_precisions (hier);






