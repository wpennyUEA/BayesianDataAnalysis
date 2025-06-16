function [mix_clean] = spm_mix_cleanup (mix,prt)
% Remove little-used clusters
% FORMAT [mix_clean] = spm_mix_cleanup (mix,prt)
%
% mix               Original model
% prt               remove clusters with pi(k) > prt
%
% mix_clean         Returned model
%

if nargin < 2, prt = 0.01; end

if mix.m == 1
    mix_clean = mix;
    return
end

mix_clean.nin = mix.nin;
mix_clean.prior = mix.prior;
mix_clean.fm = mix.fm;
mix_clean.acc = mix.acc;
mix_clean.kl_proportions = mix.kl_proportions;
mix_clean.kl_covs = mix.kl_covs;
mix_clean.kl_centres = mix.kl_centres;

pi = mix.lambda/sum(mix.lambda);
ind = find(pi>prt);
for i=1:length(ind),
    j = ind(i);
    mix_clean.state(i) = mix.state(j);
    mix_clean.lambda(i) = mix.lambda(j);
    mix_clean.gamma = mix.gamma(ind,:);
end
mix_clean.m = length(ind);

