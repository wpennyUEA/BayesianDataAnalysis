function [u,y,A] = nca_create_data (map, K, T, verbose)
% Create synthetic data for analysis
% FORMAT [u,y,A] = nca_create_data (map, K, T, verbose)
%
% map       'diff','diff2','sum','eye-x2'
% K         number of features
% T         number of data points
% verbose   default is 0
%
% u         [K x T] data to be classified
% y         [1 x T] class labels
% A         [P x K] true feature matrix

if nargin < 4 | isempty(verbose)
    verbose=0;
end

u_max = 10;
u = rand(K,T)*u_max;

switch map,
    
    case 'diff',
        A = zeros(1,K);
        A (1,1:2) = [1 -1];
        x = A*u;
        y = abs(x)<2.5;
        
    case 'diff2',
        A = zeros(2,K);
        A (1,1:2) = [-1 1];
        A (2,3:4) = [1 -1];
        x = A*u;
        y = abs(x(1,:))<4.5 & abs(x(2,:))<4.5;
        
    case 'sum',
        A = zeros(1,K);
        A(1,1:2) = [1 1];
        x = A*u;
        y = x > 10;
    
    case 'eye-x2',
        A = eye(K);
        x = A*u;
        %y = x > 10;
        sx=zeros(1,T);
        for k=1:K,
            sx=sx+x(k,:).^2;
        end
        sx=sx/(K*u_max^2);
        msx = median(sx);
        y = sx > msx;
end

if verbose
    p1=length(find(y==1))/T;
    disp(' ');
    disp(sprintf('Prob y=1 is %1.2f',p1));
end
disp(' ');