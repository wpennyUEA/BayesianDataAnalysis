
clear all
close all

% This script compares Savage-Dickey and Default JZS Baues Factors in a
% one-sample t-test when the prior is properly scaled to the data (not an arbitrary scale)

% Number of parameters (intercept only)
K=1;

N=16; % Sample size
glm.X=ones(N,K);    % One-sample design matrix

% Prior hyperparameters
glm.c0=2;   % Shape of Gamma prior
glm.b0=0.5; % Scale of Gamma prior
glm.B0=1*eye(K);    % Prior covarience
glm.w0=zeros(K,1);  % Prior mean

% Number of repetitions/simulations
R=1000;

% Preallocate
lambda=spm_gamrnd(glm.c0,glm.b0,R,1);   % Sample precisions
iB0=inv(glm.B0);    % Inverse prior covariance

for i=1:R,
    % Sample w from prior, scaled by precision
    C0 = iB0/lambda(i);
    w(:,i) = spm_normrnd(glm.w0,C0,1);
    %w(:,i) = 0;
    % Effect size
    d(i) = w(:,i)*sqrt(lambda(i)); 
    % Simulate data
    y = glm.X*w(:,i)+sqrt(1/lambda(i))*randn(N,1);
    % Compute Bayes Factors and t-statistic
    [logbf(i),logbf_Default(i),t(i)] = bayes_glm_ttest1 (y);
end

% Plot comparison of Log Bayes Factors
figure
subplot(2,2,1);
plot(logbf_Default,logbf,'x');
xlabel('LogBF Default');
ylabel('LogBF SDT');
grid on
[tmp,ind]=sort(logbf_Default);
hold on
plot(tmp,tmp,'r-');

% Plot effect size vs LogBF
subplot(2,2,2);
plot(d,logbf_Default,'bx');
hold on
grid on
plot(d,logbf,'rx');
xlabel('True Effect Size');
ylabel('LogBF');
legend('Default','SDT');

% Plot t-statisitc vs LogBF
subplot(2,2,4);
plot(t,logbf_Default,'bx');
hold on
grid on
plot(t,logbf,'rx');
xlabel('t-statistic');
ylabel('LogBF');
legend('Default','SDT');


