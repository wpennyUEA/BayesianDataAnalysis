
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

rx=5;ry=6;
S.xmin=-rx;S.xmax=rx;S.dx=0.1;
S.ymin=-ry;S.ymax=ry;S.dy=0.1;
h2 = vbcca_marginal_2D (cca,S,X1,X2);

% c2 = [-5:0.1:5];
% C1p = [-10:0.1:10];
% con.Gamma1 = [1 0]; % first dimension of x1
% con.Gamma2 = [1 0]; % first dimension of x2
% N = length(c2);
% for n=1:N,
%     [c1hat(n),gamma(:,n),p1(n,:)] = vbcca_cond_subspace (cca,c2(n),con, C1p);
% end
% figure;
% imagesc(c2,C1p,p1');
% axis xy
% xlabel('x2[1]');
% ylabel('p(x1[1]|x2[1])');
% hold on

% Predict X1[1] from both X2 variables
con.Gamma1 = [1 0]; % first dimension of x1
con.Gamma2 = eye(2); % both dimensions of x2
for n=1:size(X2,2),
    c1_both(n) = vbcca_cond_subspace (cca,X2(:,n),con);
end

% Predict X1[1] from both X2[1]
con.Gamma1 = [1 0]; % first dimension of x1
con.Gamma2 = [1 0]; 
for n=1:size(X2,2),
    c1_first(n) = vbcca_cond_subspace (cca,X2(1,n),con);
end

% Predict X1[1] from X2 [2]
con.Gamma1 = [1 0]; % first dimension of x1
con.Gamma2 = [0 1]; 
for n=1:size(X2,2),
    c1_second(n) = vbcca_cond_subspace (cca,X2(2,n),con);
end

figure
plot(X2(1,:),con.Gamma1*X1,'bo');
hold on
plot(X2(1,:),c1_both,'r.');
plot(X2(1,:),c1_first,'k.');
plot(X2(1,:),c1_second,'g.');
xlabel('X2[1]');
ylabel('X1[1]');
legend({'Original Data','Both X2','X2[1]','X2[2]'});



figure
plot(con.Gamma1*X1,c1_both,'ro');
hold on
plot(con.Gamma1*X1,c1_first,'ko');
plot(con.Gamma1*X1,c1_second,'go');
xlabel('Data');
ylabel('Fit');
legend({'Both X2','X2[1]','X2[2]'});


% % Plot original data, but adjusted for contrasts
% [X1adj,X2adj,c1adj] = vbcca_adjusted (cca,con,X1,X2);
% plot(X2adj,X1adj,'wo');
% %plot(con.Gamma2*X2,con.Gamma1*X1,'wo');
% 
% 
% 
% figure;
% plot(c2,c1hat,'r');
% hold on
% plot(X2adj,c1adj,'bo');
% legend({'Pred from X2(1)','Pred from X2(2)'});
% 
