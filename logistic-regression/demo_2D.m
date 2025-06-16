
clear all
close all

disp('Demo of Bayesian Logistic Regression with 2D inputs');

N=120; % Number of data points

opt.verbose=1;

% Generate data from mixture model
mix.m=2;
mix.state(1).m=[1,1];
mix.state(2).m=[3,3];
for i=1:2,
    mix.state(i).C=eye(2);
    mix.state(i).prior=1/2;
end

[x,label] = spm_samp_mix (mix,N);

if opt.verbose
    figure
    col={'rx','bx'};
    for i=1:2,
        ind=find(label==i);
        plot(x(ind,1),x(ind,2),col{i});
        hold on
    end
end

% Randomly permute data
ind = randperm(N);
x = x(ind,:);
label = label(ind,:);

opt.diag=0;
opt.verbose=1;

t = label-1;
M = blr_fit (x,t,opt);

blr_plot_boundary (M,x);
