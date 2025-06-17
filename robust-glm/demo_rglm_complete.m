% This script compares parameter estimation errors between a standard GLM using ordinary least squares, and then a
% robust GLM with 1 and 2 mixture components

% Generate data from ARMOG model
N=351;
x1=spm_boxcars(N,1,10); % Design regressor
x2=ones(N,1);   % Intercept
X=[x1,x2]; 
beta=[1 1]';    %True parameters

% Define Gaussian-Mixture noise parameters (2 components)
m=2;
mix.m=m;
mix.state(1).prior=0.73;
mix.state(2).prior=0.27;
mix.state(1).m=0;
mix.state(2).m=0;
% Variance of the components
mix.state(1).C=2.4^2;
mix.state(2).C=8.4^2;

% Sample noise from Gaussian mixture
[noise,gamma_true]=spm_samp_mix(mix,N);
new_index=randperm(N);
noise=noise(new_index);
gamma_true=gamma_true(new_index);

% Generate data with mixture noise
y=X*beta + noise;

% Fit robust GLM with 1 and 2 mixture components
m=1;
rglm1 = spm_rglm (y,X,m);
m=2;
rglm2 = spm_rglm (y,X,m);
% Fit standard GLM (ordinary least squares)
w_ml=pinv(X)*y;

% Plot true signal and observed data
figure
plot(X*beta);
hold on;
plot(y,'r');
legend('Signal','Data');

% Plot histogram of noise
figure
hist(noise,20);
title('Error histogram');

% Plot estimated vs true class labels for mixture components
figure
plot(rglm2.posts.gamma(1,:),'r');
hold
plot(gamma_true);
axis([0 N -0.1 2.1]);
title('Class labels');
legend('Estimated','True');

% Calculate error norms for parameter estimates
d_ml=norm(w_ml-beta);
d_vb=norm(rglm2.posts.w_mean-beta);

% Model averaging (based on model evidence)
f=[rglm1.fm,rglm2.fm];
f=f-mean(f);    % Normalise
ef=exp(f);
pm=ef./sum(ef); % Posterior probability of each model
w_ma=pm(1)*rglm1.posts.w_mean + pm(2)*rglm2.posts.w_mean;
d_ma=norm(w_ma-beta);

% Display results
disp(' ');
disp(sprintf('Error for GLM=%1.3f',d_ml));
disp(sprintf('Error for RGLM=%1.3f',d_vb));

disp(' ');
err_ratio=d_ml/d_vb;
disp(sprintf('Ratio MSE_GLM/MSE_RGLM=%1.3f',d_ml/d_vb));

% Calculate negative predictive log-likelihoods (cost)

% 1. ML Gaussian likelihood
e_ml=y-X*w_ml;  
s_ml=std(e_ml); 
p_ml=spm_Npdf(e_ml,0,s_ml^2);   % Evaluate Gaussian PDF at residuals
E=-log(p_ml);   % Negative log-likelihood
% Compute theoretical NLL curve for a range of error values
el=[-4*s_ml:0.1:4*s_ml];
pl=spm_Npdf(el,0,s_ml^2);
El=-log(pl);

% 2. Variational Bayes (Gaussian Mixture)
p_vb=zeros(N,1);
for s=1:rglm2.m,
    p_g=spm_Npdf(e_ml,0,rglm2.posts.variances(s)); % Evaluate Gaussian PDF at residuals with component variance
    p_vb=p_vb+rglm2.posts.pi(s)*p_g;    % Weight by the component's posterior probability
end
E_vb=-log(p_vb);    % Negative log-likelihood
% Compute theoretical NLL curve for a range of error values
p_v=zeros(size(el));
for s=1:rglm2.m,
    p_g=spm_Npdf(el,0,rglm2.posts.variances(s));
    p_v=p_v+rglm2.posts.pi(s)*p_g;
end
Ev=-log(p_v);

% Plot theoretical NLL
figure
plot(el,El,'LineWidth',2);
hold on
plot(el,Ev,'r','LineWidth',2);
[ee,ei]=sort(e_ml);
% Plot NLL values
plot(e_ml,E,'.','MarkerSize',20);
plot(e_ml,E_vb,'r.','MarkerSize',20);

set(gca,'FontSize',18);
legend('Gaussian','MoG');
xlabel('Error, e_n');
ylabel('-ln p(e_n)');
title('Cost');

% Compute Z-values (effect size over uncertainty)
z1=rglm1.posts.w_mean(1)/sqrt(rglm1.posts.w_cov(1,1));
z2=rglm2.posts.w_mean(1)/sqrt(rglm2.posts.w_cov(1,1));

disp(sprintf('Z value for GLM = %1.2f',z1));
disp(sprintf('Z value for RGLM = %1.2f',z2));


