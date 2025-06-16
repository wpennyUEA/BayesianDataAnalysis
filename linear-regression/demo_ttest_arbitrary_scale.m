
clear all
close all

disp('One sample t-test demo - data is on arbitrary scale ...');
disp('Compare Savage Dickey with T priors against Default approach');
disp(' ');

% Agreement gets better with increasing N
N=16;
X=ones(N,1);

sd=6;
m=[-6:0.1:6];

for i=1:length(m),
    d(i) = m(i)/sd;  % Effect size
    y = X*m(i)+sd*randn(N,1);
    [logbf(i),logbf_jzs(i),t(i)] = bayes_glm_ttest1 (y);
end

figure
subplot(2,2,1);
plot(logbf_jzs,logbf,'x');
xlabel('LogBF Default');
ylabel('LogBF SDT');
grid on
[tmp,ind]=sort(logbf_jzs);
hold on
plot(tmp,tmp,'r-');

subplot(2,2,2);
plot(d,logbf_jzs,'bx');
hold on
grid on
plot(d,logbf,'rx');
xlabel('True Effect Size');
ylabel('LogBF');
legend('Default','SDT');

subplot(2,2,4);
plot(t,logbf_jzs,'bx');
hold on
grid on
plot(t,logbf,'rx');
xlabel('t-statistic');
ylabel('LogBF');
legend('Default','SDT');


