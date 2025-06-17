
clear all
close all

% This script compares the Savage-Dickey test using a t-prior against the
% Default JZS Bayes Factor approach in one-sample t-tests (data is on an arbitrary scale)

N=16;   % Sample size - agreement gets better with increasing N
X=ones(N,1);    % Design matrix: single regressor (intercept only)
sd=6;   % Noise SD
m=[-6:0.1:6];   % Range of true means to simulate


for i=1:length(m),
    % Effect size
    d(i) = m(i)/sd; 
    % Simulate data
    y = X*m(i)+sd*randn(N,1);
    % Compute Bayes Factors and t-statistic
    [logbf(i),logbf_jzs(i),t(i)] = bayes_glm_ttest1 (y);
end

% Plot SDT vs Default Bayes Factor
figure
subplot(2,2,1);
plot(logbf_jzs,logbf,'x');
xlabel('LogBF Default');
ylabel('LogBF SDT');
grid on
[tmp,ind]=sort(logbf_jzs);
hold on
plot(tmp,tmp,'r-');

% Plot effect size vs LogBF
subplot(2,2,2);
plot(d,logbf_jzs,'bx');
hold on
grid on
plot(d,logbf,'rx');
xlabel('True Effect Size');
ylabel('LogBF');
legend('Default','SDT');

% Plot t-statistic vs LogBF
subplot(2,2,4);
plot(t,logbf_jzs,'bx');
hold on
grid on
plot(t,logbf,'rx');
xlabel('t-statistic');
ylabel('LogBF');
legend('Default','SDT');


