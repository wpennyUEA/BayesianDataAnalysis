function [R] = robust_reg (y,x,verbose)
% Robust univariate regression
% FORMAT [R]  = robust_reg (y,x,verbose)
%
% y         [N x 1] dependent variable
% x         [N x 1] independent variable
% verbose   text and figure outputs (default = 1)
%
% R     .y  dependent variable with outliers cleaned up
%       .r  correlation coefficient for "cleaned" data

try verbose = verbose; catch verbose = 1; end

y = y(:); x = x(:);
N = length(x);
Ncheck = length(y);
if ~(N==Ncheck)
    disp('Error in robust_unireg: vectors not same length');
    return
end

X = [x,ones(N,1)];
rglm1 = spm_rglm(y(:),X,1);
[rglm2,yclean] = spm_rglm(y(:),X,2);
if rglm1.fm > rglm2.fm
    if verbose
        disp(sprintf('Single-Component Gaussian Error model preferred, LogBF=%1.2f',rglm1.fm-rglm2.fm));
    end
    R.r = corr(x,y);
    R.yclean = y;
else
    if verbose
        disp(sprintf('Two-Component Gaussian Error model preferred, LogBF=%1.2f',rglm2.fm-rglm1.fm));
    end
    
    for m=1:2,
        disp(' ');
        disp(sprintf('Error component %d:',m));
        disp(sprintf('Proportion of samples = %1.2f',rglm2.posts.pi(m)));
        disp(sprintf('Error SD = %1.2f',sqrt(rglm2.posts.variances(m))));
    end
    R.r = corr(x,yclean);
end



if verbose
    figure
    plot(x,y,'x');
    hold on
    beta = pinv(X)*y;
    Xr = [min(x) 1;max(x) 1];
    yhat1 = Xr*beta;
    plot([min(x),max(x)],[yhat1(1) yhat1(2)],'r-');
    
    if rglm2.fm > rglm1.fm
        yhat2 = Xr*rglm2.posts.w_mean;
        plot([min(x),max(x)],[yhat2(1) yhat2(2)],'k-');
        legend('Data','StandardFit','RobustFit');
    
        title(sprintf('Robust correlation, r = %1.2f', R.r));
    else
        legend('Data','StandardFit');
        title(sprintf('Correlation, r = %1.2f', R.r));
    end
end



