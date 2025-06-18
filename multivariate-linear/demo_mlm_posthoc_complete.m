
clear all
close all

% This script demonstrates posthoc model comparison for spm_mlm_bayes.

% Set dimensions
N=100;
d=10;
p=5;

% Generate design matrix 
x=randn(N,p);
% Generate true regression coefficients. In this example only inputs 1,2 and 3 are predictive of outputs
W=randn(p,d);
W(4:5,:)=0;
% Generate noise
e=2*randn(N,d);

% Generate observed data (linear model + noise)
y=x*W+e;

% Setup options for spm_mlm_bayes
options.pr = 'input';   % Input prior shrinks connections from irrelevant inputs to zero
options.verbose = 1;

% Run Bayesian MLM estimation
evalc('mlm = spm_mlm_bayes (y,x,options);');

% Display Bayesian regression coefficients (posterior means)
figure
imagesc(mlm.wmean);
colormap gray
colorbar
ylabel('Inputs');
xlabel('Outputs');
title('Bayes Regression Coefficients');

% Display Maximum Likelihood regression coefficents for comparison
figure
imagesc(mlm.wml);
colorbar
ylabel('Inputs');
xlabel('Outputs');
colormap(gray);
title('ML Regression Coefficients');

% Post-hoc inference based on Savage-Dickey tests
disp(' ');
disp('Posthoc questions:');
disp(' ');

disp('Is regression coefficient 1,1 non zero ?');
con=zeros(p,d);
con(1,1)=1; % Define contrast vector selecting coefficient (1,1)
con_vec=con(:)';
disp('Log Evidence in favour of this hypothesis:');
logbf = spm_mlm_posthoc (mlm,con_vec)

disp('Is regression coefficient 4,6 non zero ?');
con=zeros(p,d);
con(4,6)=1;
con_vec=con(:)';
disp('Log Evidence in favour of this hypothesis:');
logbf = spm_mlm_posthoc (mlm,con_vec)

disp('Are regression coefficients 1,1 2,1 and 3,1 non zero ?');
w=zeros(p,d);
w(1,1)=1;w(2,1)=1;w(3,1)=1;
con_vec1 = spm_mlm_makecon (mlm,w);
disp('Log Evidence in favour of this hypothesis:');
logbf = spm_mlm_posthoc (mlm,con_vec1)

disp('If we define w1 = regression coefficients 1,1 2,1 and 3,1');
disp('and w2 = regression coefficients 1,2 2,2 and 3,2');
disp('Is w1 different to w2 ?');
w=zeros(p,d);
w(1,2)=1;w(2,2)=1;w(3,2)=1;
con_vec2 = spm_mlm_makecon (mlm,w);
con_vec_diff=con_vec2-con_vec1; % Contrast difference between w2 and w1
disp('Log Evidence in favour of this hypothesis:');
logbf = spm_mlm_posthoc (mlm,con_vec_diff)
