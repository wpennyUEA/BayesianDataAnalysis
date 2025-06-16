function [mix,px,newcause] = igmm_update (mix,x)
% Update mixture model given new data point
% FORMAT [mix,px,newcause] = igmm_update (mix,x)
%
% mix       mixture model data structure
% x         [D x 1] data vector
%
% mix       updated model
% px        probability of data point
% newcause  1 for new cause created, 0 otherwise
%
% [1] Pinto & Engel (2015) A Fast Incremental Gaussian Mixture 
% Model. PLoS One.

for j=1:mix.M,
    err(:,j) = x-mix.state(j).m;
    d2(j) = err(:,j)'*mix.state(j).Lambda*err(:,j);
end

like0=(2*pi)^(-mix.D/2);
% Update Components
for j=1:mix.M
    %like1 =dbquit
    like(j) = like0*(mix.state(j).detC)^(-0.5)*exp(-0.5*d2(j));
    post(j) = mix.prior(j)*like(j);
end
px=sum(post);
post=post/px;

newcause = 0;
if min(d2) > mix.chi2
    mix = igmm_create(mix,x);
    newcause = 1;
    return
end

mix.v = mix.v+1;
mix.sp = mix.sp+post;
mix.prior = mix.sp/sum(mix.sp); 
w=post./mix.sp;

debug=0;

for j=1:mix.M,
    dmu = w(j)*err(:,j);
    mix.state(j).m = mix.state(j).m + dmu;
    new_err = x-mix.state(j).m;
    
    % Equation 20 in [1]
    Lambda = mix.state(j).Lambda;
    w1 = 1-w(j);
    term1 = Lambda/w1;
    num = (w(j)/(w1^2)) * Lambda *new_err*new_err'*Lambda;
    denom1 = 1 + (w(j)/w1) * new_err'*Lambda*new_err;
    Lambda_bar = term1 - num/denom1;
    
    % Equation 21 in [1]
    num = Lambda_bar*dmu*dmu'*Lambda_bar;
    denom2 = 1 - dmu'*Lambda_bar*dmu;
    mix.state(j).Lambda = Lambda_bar + num/denom2;
    
    if debug
        % Check implementation using equation 11  
        % for covariance update 
        C = inv(Lambda);
        Cnew = w1*C + w(j)*new_err*new_err' - dmu*dmu';
        Lambda_new = inv(Cnew);
        %disp(mix.state(j).Lambda);
        %disp(Lambda_new);
        Cerr=sum(sum(Lambda_new-mix.state(j).Lambda));
        if Cerr > 1e-6 
            disp('Error in Lambda update');
            keyboard
        end
    end
    
    % Equations 25 and 26 in [1]
    term1=(1-w(j))^mix.D;
    Cbar = term1*mix.state(j).detC*denom1;
    mix.state(j).detC = Cbar*denom2;
    
    if debug
        Cerr=mix.state(j).detC-det(Cnew);
        if Cerr > 1e-6 
            disp('Error in Determinant update');
            keyboard
        end
    end
    if mix.state(j).detC < 0 
        disp('Error in igmm_update');
        disp('Increase s0 parameter');
        return
    end
end


