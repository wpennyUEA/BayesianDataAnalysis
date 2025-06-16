
clear all
close all

disp('This script compares e.g. Savage Dickey Bayes Factors');
disp('to Default Bayes Factors');

N=100;  % Number of data points
P=4;    % Number of regressors

% Number of simulations
Reps=50;

% Savage-Dickey/MultiFit Bayes Factors based on T or Gaussian distributions ?
method='T';
%method='G';

% Average proportion of variance explained 
% for simulations where alt is true
R2_avg = 0.5;

% Which two types of LogBFs to compare ?
logBF_type={'Savage-Dickey','MultiFit','Default'};
ix=3;
iy=1;

general_design=1;
if general_design
    % Design matrix with correlated non-zero mean regressors
    rx=0.7; mx=3;
    X(:,1)=randn(N,1)+mx;
    for i=2:P-1,
        X(:,i)=X(:,i-1)+rx*randn(N,1)+mx;
    end
    X(:,P)=ones(N,1);
else
    % Design matrix is from a DCT which has orthogonal, zero-mean columns
    X=spm_dctmtx(N,P);
    % Final column of design matrix must be 1's
    X(:,1)=X(:,P);
    X(:,P)=ones(N,1);
end

glm.X=X;
glm.pE=zeros(P,1);
glm.pC=eye(P);

% Work out appropriate observation noise level to get R2_avg
w_true = spm_normrnd(glm.pE,glm.pC,Reps);
y_true = X*w_true;
vy = mean(std(y_true).^2);
ve = vy*(1-R2_avg)/R2_avg;
glm.Ce = ve*eye(N);

% Run for (h=1) null is true and (h=2) alt is true
alt_true=[0,1];
hname={'Null True','Alt True'};

for h=1:2,
    for r=1:Reps,
        if alt_true(h), wr=w_true(:,r); else wr = zeros(P,1); end
        e = sqrt(ve)*randn(N,1);
        y = X*wr+e;
        
        [logBF,glm_post] = bayes_glm_regression (X,y,method);
        logbf_x(h,r)= logBF(ix);
        logbf_y(h,r) = logBF(iy);
        R2(h,r) = glm_post.R2;
    end
    subplot(2,2,h);
    plot(logbf_x(h,:),logbf_y(h,:),'x');
    hold on
    xlabel(logBF_type{ix});
    ylabel(logBF_type{iy});
    grid on
    [tmp,ind]=sort(logbf_x(h,:));
    hold on
    plot(tmp,tmp,'r-');
    title(hname{h});
    
    subplot(2,2,h+2);
    plot(R2(h,:),logbf_x(h,:),'rx');
    hold on
    grid on
    plot(R2(h,:),logbf_y(h,:),'x');
    legend({logBF_type{ix},logBF_type{iy}});
    xlabel('Proportion Variance Explained');
    ylabel('LogBF');
end

