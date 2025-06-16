
clear all
close all

disp('Demo of Bayesian Logistic Regression with Ionosphere data');

load iono_data.mat

t = iono_data(:,35) - 1;
x = iono_data(:,[1,3:34]);  % 2nd variable has zero SD, so ignore
[N,P] = size(x);

opt.verbose=1;
M = blr_fit (x,t,opt);

figure
plot(M.z(1:end-1));
xlabel('Parameter');
ylabel('Z-score');
grid on
hold on
plot([1 P],[1.96,1.96],'r-');
plot([1 P],[-1.96,-1.96],'r-');

% Compare models
model(1).x=[];
model(2).x=x;
model(3).x=x(:,26);
ind=find(abs(M.z(1:P))>1.96);
model(4).x=x(:,ind);
name={'const','All','feat26','All-Sig'};

F = blr_compare(model,name,t,opt);
