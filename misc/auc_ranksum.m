function [auc,pval] = auc_ranksum (x,y)
% Compute AUC and pvalue using ranksum
% FORMAT [auc,pval] = auc_ranksum (x,y)
%
% x     x scores
% y     y scores

if mean(x) > mean(y)
    % Swap x and y
    tmp = y;
    y = x;
    x = tmp;
end
    
Nx = length(x);
Ny = length(y);
[pval,h,stats] = ranksum(y,x);

% Mann-Whitney U 
U = stats.ranksum-0.5*Ny*(Ny+1);
auc = U/(Ny*Nx);