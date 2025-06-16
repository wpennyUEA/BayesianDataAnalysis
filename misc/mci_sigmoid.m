function [y] = mci_sigmoid (x)
% Sigmoid function
% FORMAT [y] = mci_sigmoid (x)

y = 1./(1+exp(-x));

end