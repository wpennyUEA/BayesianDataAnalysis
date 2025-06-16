function [y,a,ymod,sigma] = bsr_output (bsr,U)
% Return output of softmax regression model
% FORMAT [y,a,ymod,sigma] = bsr_output (bsr,U)
%
% bsr       data structure - see bsr_fit.m
% U         inputs
%
% y         outputs
% a         activations
% ymod      moderated outputs
% sigma     standard deviation of activations

% Restrict -alim < a < alim
alim=100;

class=bsr.class;
K=length(class);
for k=1:K,
    a(:,k)=U*class(k).w;
end
a=min(a,alim);
a=max(a,-alim);

ea=exp(a)+eps;
sea=sum(ea,2)+eps;
y = ea./(sea*ones(1,K));

if nargout > 2
    
    % The below implements eq 9.61 in Nabney 2003
    for k=1:K,
        try iR=class(k).iR; catch iR=inv(class(k).R); end
        
        % Is diag(U*iR*U') quicker?
        N=size(U,1);
        for n=1:N,
            s2(n,k)=U(n,:)*iR*U(n,:)';
        end
        h(:,k)=hfun(s2(:,k));
        amod(:,k)=h(:,k).*a(:,k);
    end
    amod=min(amod,alim);
    amod=max(amod,-alim);

    eamod=exp(amod)+eps;
    seamod=sum(eamod,2)+eps;
    ymod = eamod./(seamod*ones(1,K));
    sigma = sqrt(s2);
end

end

%------------------------------------------------------------------
function [h] = hfun (s2)

h0=1+(s2*pi)/8;
h=h0.^(-0.5);

end