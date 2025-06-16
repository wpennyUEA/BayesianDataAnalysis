
clear all
%close all

% Number of data points
N=100;

% Number of regressors
p=5;

r2=[0:0.02:0.8];
for i=1:length(r2),
    bf10 = linregbf(r2(i),N,p);
    logbf_LR(i) = log(bf10);
    if p==1
        bf10 = corrbf(sqrt(r2(i)),N);
        logbf_Corr(i) = log(bf10);
    end
end

figure
plot(r2,logbf_LR,'b-');
grid on
xlabel('Prop Variance Explained, r2');
ylabel('LogBF alt vs null');
if p==1
    hold on
    plot(r2,logbf_LR,'r-');
    legend('LR','Corr');
end
