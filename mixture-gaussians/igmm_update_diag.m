function [mix] = igmm_update_diag (mix,x)
% Update diagonal cov mixture model given new data point
% FORMAT [mix] = igmm_update_diag (mix,x)
%
% mix       mixture model data structure
% x         [D x 1] data vector
%
% Pinto & Engel (2015) A Fast Incremental Gaussian Mixture 
% Model. PLoS One.

for j=1:mix.M,
    err(:,j) = x-mix.state(j).m;
    d2(j) = err(:,j)'*mix.state(j).Lambda*err(:,j);
end

if min(d2) > mix.chi2
    mix = igmm_create(mix,x);
    return
end

like0=(2*pi)^(-mix.D/2);

% Update Components
for j=1:mix.M
    like(j) = like0*(mix.state(j).detC)^(-0.5)*exp(-0.5*d2(j));
    post(j) = mix.prior(j)*like(j);
end
post=post/sum(post);

mix.v = mix.v+1;
mix.sp = mix.sp+post;
mix.prior = mix.sp/sum(mix.sp); 
w=post./mix.sp;

debug=1;

for j=1:mix.M,
    dmu = w(j)*err(:,j);
    new_err = x-mix.state(j).m;
    
    C = diag(1./diag(mix.state(j).Lambda));
    % Equation 11 in [1]
    Cnew = (1-w(j))*C + w(j)*new_err*new_err' - dmu*dmu';
    
    mix.state(j).m = mix.state(j).m + dmu;
    mix.state(j).Lambda = diag(1./diag(Cnew));
    mix.state(j).detC = prod(diag(Cnew));
end


