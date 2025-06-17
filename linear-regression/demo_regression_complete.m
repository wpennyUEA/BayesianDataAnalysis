
clear all
close all

% This script compares different Bayes Factors from different methods in
% Bayesian regression: Savage Dickey, Multifit and Default Bayes Factors

N=100;  % Number of data points
P=4;    % Number of regressors

% Number of simulations
Reps=50;

% Method for computing posteriors: T or Gaussian distributions?
method='T';
%method='G';

% Average proportion of variance explained 
% for simulations where alt is true
R2_avg = 0.5;

% Which two Bayes Factor types to compare? 1=Savage-Dickey, 2=MultiFit,
% 3=Default
logBF_type={'Savage-Dickey','MultiFit','Default'};
ix=3;
iy=1;

% Create design matrix (either general or orthogonal)
general_design=1;
if general_design
    % Create correlated regressors with non-zero mean (general design)
    rx=0.7; mx=3;
    X(:,1)=randn(N,1)+mx;
    for i=2:P-1,
        X(:,i)=X(:,i-1)+rx*randn(N,1)+mx;
    end
    X(:,P)=ones(N,1);   % Intercept
else
    % Use DCT regressors with zero-mean (orthogonal design))
    X=spm_dctmtx(N,P);
    X(:,1)=X(:,P);  % Final column becomes first
    X(:,P)=ones(N,1);   % Intercept
end

% Set up GLM prior parameters
glm.X=X;
glm.pE=zeros(P,1);
glm.pC=eye(P);

% Work out appropriate observation noise level to get R2_avg
w_true = spm_normrnd(glm.pE,glm.pC,Reps);   % True weights from the prior
y_true = X*w_true;
vy = mean(std(y_true).^2);  % Estimate variance in true signal across reps
ve = vy*(1-R2_avg)/R2_avg;  % Set noise variance
glm.Ce = ve*eye(N); % Error covariance

% Run simulations: (h=1) null is true and (h=2) alt is true
alt_true=[0,1];
hname={'Null True','Alt True'};
for h=1:2,
    for r=1:Reps,
        if alt_true(h), wr=w_true(:,r); % True weights
        else wr = zeros(P,1);   % All weights zero (null)
        end
        e = sqrt(ve)*randn(N,1);  % Generate noise
        y = X*wr+e;  % Simulated data

        % Run Bayesian GLM and store Bayes factors
        [logBF,glm_post] = bayes_glm_regression (X,y,method);
        logbf_x(h,r)= logBF(ix);
        logbf_y(h,r) = logBF(iy);
        R2(h,r) = glm_post.R2;
    end

    % Plot comparison of selected LofBF types
    subplot(2,2,h);
    plot(logbf_x(h,:),logbf_y(h,:),'x');
    hold on
    xlabel(logBF_type{ix});
    ylabel(logBF_type{iy});
    grid on
    % Add identity line to see agreement
    [tmp,ind]=sort(logbf_x(h,:));
    hold on
    plot(tmp,tmp,'r-');
    title(hname{h});
    
    % Plot relationship between R2 and LogBFs
    subplot(2,2,h+2);
    plot(R2(h,:),logbf_x(h,:),'rx');
    hold on
    grid on
    plot(R2(h,:),logbf_y(h,:),'x');
    legend({logBF_type{ix},logBF_type{iy}});
    xlabel('Proportion Variance Explained');
    ylabel('LogBF');
end

