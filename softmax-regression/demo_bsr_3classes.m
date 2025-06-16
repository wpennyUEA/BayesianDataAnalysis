
clear all
close all

disp('Three class, 2D inputs demo of Bayesian Softmax Regression');

N=120; % Number of data points
Dnoise=0; % Number of additional noisy inputs (distractors)
disp(sprintf('Number of spurious predictors added = %d',Dnoise));

opt.verbose=1;

% Generate data from mixture model
mix.m=3;
mix.state(1).m=[1,1]';
mix.state(2).m=[1,5]';
mix.state(3).m=[3,3]';
for i=1:3,
    mix.state(i).C=eye(2);
    mix.state(i).prior=1/3;
end

[x,label] = spm_samp_mix (mix,N);

if opt.verbose
    figure
    col={'rx','bx','kx'};
    for i=1:3,
        ind=find(label==i);
        plot(x(ind,1),x(ind,2),col{i});
        hold on
    end
end

% Randomly permute data
ind = randperm(N);
x = x(ind,:);
label = label(ind,:);

% Add spurious random predictors
x = [x,1.65*randn(N,Dnoise),ones(N,1)];

opt.diag=0;
opt.verbose=1;
tic; bsr = bsr_fit (x,label,opt); toc

[y,a] = bsr_output (bsr,x);

[logbf,z] = bsr_savage_dickey(bsr);
disp(' ');
disp('Log BF, row is class(k), column is feature (d):');
disp(logbf)
disp(' ');
disp('z, row is class(k), column is feature (d):');
disp(z)
disp(' ');
disp('Sum LogBF over classes:');
disp(sum(logbf,1));

[y,a,ymod] = bsr_output (bsr,x);

if opt.verbose
    figure
    col={'rx','bx','kx'};
    for i=1:3,
        [tmp,assign]=max(y');
        ind=find(assign'==i);
        plot(x(ind,1),x(ind,2),col{i});
        hold on
    end
    if Dnoise==0
        % If Dnoise>0 there will be > 2 inputs
        bsr_plot_boundary (bsr,x);
    end
       
end