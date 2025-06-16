function [sens,spec,stats] = sensitivity_pcut (x,y,pcut)
% Compute sensitivity and specificity at CDF cutoff
% FORMAT [sens,spec,stats] = sensitivity_pcut (x,y,pcut)
%
% x     [1 x N1] data from category 1
% y     [1 x N2] data from category 2
% pcut  cutoff percentile e.g. 0.95
%
% sens  sensitivity
% spec  specificity
% stats .tn,.tp,.fn,.fp true/false positive negatives
%       .xpos,.ypos indices of above threshold data

allh = [x,y];
[tmp,ind]=sort(allh);
Nh=length(allh);
Ncut=round(pcut*Nh);
hcut=tmp(Ncut);

stats.xpos = find(x>hcut);
stats.ypos = find(y>hcut);

stats.tp = length(stats.xpos);
stats.tn = length(find(y<hcut));
stats.fp = length(stats.ypos);
stats.fn = length(find(x<hcut));
sens = stats.tp/(stats.tp+stats.fp);
spec = stats.tn/(stats.tn+stats.fn);