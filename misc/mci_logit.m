function [x] = mci_logit (y,a)
% Logit Function
% FORMAT [x] = mci_logit (y,a)

x = log(y./(a-y));

end