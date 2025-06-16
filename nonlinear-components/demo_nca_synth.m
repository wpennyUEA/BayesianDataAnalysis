
clear all
close all

%map='eye-x2';
%map='diff2';
map='diff';
%map='sum';

K = 4; % Number of input features
T = 200; % Number of data points

% Which NCA algorithm ?
algorithm='BP';

[u,y,A] = nca_create_data (map, K, T);

M.verbose = 1;
switch algorithm
    case 'ML',
        disp('Maximum Likelihood Algorithm');
        %M.p = 1;
        %M.opt = 'FixedStepSize';
        M = nca_org (u', y', M);
        
    case 'BR',
        disp('Bayesian Regularization Algorithm');
        M.prune_its_min=1000;
        M = nca_prune (u', y', M);
        
    case 'BP',
        disp('Bayesian Pruning Algorithm');
        M = nca_prune (u', y', M);
        
    otherwise
        disp('Error in demo_nca_synth.m: unknown algorithm');
        return
        
end

figure
subplot(2,1,1);
hinton(A);
title('True Feature Matrix');
subplot(2,1,2);
hinton(M.A);
title('Estimated Feature Matrix');
