
clear all
close all

disp('Create two 2D data sources with single latent variable and two clusters');
disp('Estimate factor matrices using MCCA');

% Number of data points
N = 100;

% Latent variable
z = randn(1,N);

% Factor matrices  
W1{1} = [0.5 -0.3]'; % for first
W2{1} = [1 -2]';
W1{2} = [-0.3 0.5]'; % and second cluster
W2{2} = [-2 1]';

% Means 
mu1{1} = [-3 3]'; % for first  
mu2{1} = [-2 2]';
mu1{2} = [2 -2]'; % and second cluster
mu2{2} = [3 -3]';

% Observation noise SD
sigma = 0.5;

% Generate data sources
X1 = []; X2 = [];
for m = 1:2,
    X1new = W1{m}*z + mu1{m}*ones(1,N) + sigma*randn(2,N);
    X2new = W2{m}*z + mu2{m}*ones(1,N) + sigma*randn(2,N);
    X1 = [X1, X1new];
    X2 = [X2, X2new];
end

% Call MCCA algorithm
cca = vbcca (X1,X2,1,2);

for m=1:2,
    disp(sprintf('CLUSTER %d',m));
    disp('Data source 1:');
    disp('True W    Estimated W');
    disp([W1{m},cca.W1{m}]);
    disp('True mean    Estimated mean');
    disp([mu1{m},cca.mu1{m}]);
    
    disp('Data source 2:');
    disp('True W    Estimated W');
    disp([W2{m},cca.W2{m}]);
    disp('True mean    Estimated mean');
    disp([mu2{m},cca.mu2{m}]);
end

figure
plot(cca.Fhist);
ylabel('Model Evidence');
xlabel('Number of iterations');
grid on

r=5;
S.xmin=-r;S.xmax=r;S.dx=0.1;
S.ymin=-r;S.ymax=r;S.dy=0.1;
h2 = vbcca_marginal_2D (cca,S,X1,X2);

plot_conditional = 1;
if plot_conditional
    % Interactive plots of conditional density, p(X1|x2) 
    % for various x2 ...
    figure
    X2 = [0 0; 1 1; 2 2; 3 3]';
    for i=1:4,
        x2 = X2(:,i);
        plot(h2,x2(1),x2(2),'rx');
        vbcca_conditional_2D (cca,S,x2);
        pause
    end
    
    figure
    X2 = [1 1; 1.25 0; 1.3 -0.5; 1.3 -0.6; 1.35 -0.7; 1.5 -1; 2 -1.5; 2.5 -2; 3 -3]';
    for i=1:size(X2,2),
        x2 = X2(:,i);
        plot(h2,x2(1),x2(2),'rx');
        vbcca_conditional_2D (cca,S,x2);
        pause
    end
end


