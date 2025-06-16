function [M,U,Xfull] = linear_model (Nobs,lambda)
% Set up data structures for linear model
% FORMAT [M,U,Xfull] = linear_model (Nobs,lambda)
%
% Nobs      number of data points
% lambda    noise precision
%
% M         model structure
% U         U.X is the design matrix

% Interval
T=100;

x0 = ones(Nobs,1);
x1 = rand(Nobs,1)*T;
x1 = sort(x1);
U.t = x1;

x1 = x1-mean(x1);
x1 = x1/std(x1);
x2 = x1.^2;
x2 = x2-mean(x2);
x2 = x2/std(x2);

U.X=[x0,x1,x2];
U.names={'Offset','Linear','Quadratic'};
            
sigma_e=sqrt(1/lambda);
M.Ce=1/lambda;

Np=size(U.X,2);
M.pE=zeros(Np,1);
M.pC=10*eye(Np);

