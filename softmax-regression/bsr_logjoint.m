function [J] = bsr_logjoint(bsr,U,c,gamma)
% Compute log joint prob of likelihood and prior
% FORMAT [J] = bsr_logjoint(bsr,U,c,gamma)
%
% bsr       data structure - see bsr_fit
% U         inputs
% c         class labels

class = bsr.class;
y = bsr_output(bsr,U);
if nargin>3
    J1 = sum(sum(gamma.*c.*log(y)));
else
    J1 = sum(sum(c.*log(y)));
end

K=length(class);
J2=0;
for k=1:K,
    dw = class(k).w-class(k).w0;
    J2 = J2 -0.5*dw'*class(k).R0*dw;
end

J = J1 + J2;