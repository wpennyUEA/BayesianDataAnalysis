function [obj] = vbcca_index_init (obj, X, Y, indices)
% Initialise using indices (assignment of data points to clusters)
% FORMAT [obj] = vbcca_index_init (obj, X, Y, indices)
%
% X         dX X N datamatrix
% Y         dY X N datamatrix
% indices   [N x 1] assignment vector
%
% output:
% obj - Model with initialised parameters.

% Default parameter settings

data = [X' Y'];
obj.N = min(cols(X), cols(Y));

%centers = zeros(obj.M, obj.dX + obj.dY);


%Initialise posterior means from centers:
% obj.eX_mu = centers(:,1:obj.dX)';
% obj.eY_mu = centers(:,(obj.dX + 1):(obj.dX + obj.dY))';

obj.r = zeros(obj.N, obj.M);

%Initialise projections & data specific variance from prob cca
%MLE results
for k=1 : obj.M

    clusterDataX = data(indices == k, 1:obj.dX);
    clusterDataY = data(indices == k,(obj.dX + 1):(obj.dX + obj.dY));

    %Hard clustering
    obj.r(indices == k, k) = 1;

    %Cluster sizes
    obj.pi(k) = sum(indices == k) / obj.N;

    % Means
    obj.eX_mu(:,k) = mean(clusterDataX);
    obj.eY_mu(:,k) = mean(clusterDataY);

    % Option for initialising with the classical CCA:
    % -------------------------------------------------
    %Data covariances
    xCov = cov(clusterDataX);
    yCov = cov(clusterDataY);

    %Normal CCA, Check Prob CCA paper for formula descriptions
    [A,B,co] =  canoncorr(clusterDataX, clusterDataY);
    ccaDim = obj.D;%min(rank(A'), rank(B'));
    codiag = diag(sqrt( co(1 : ccaDim)));

    W1 = [xCov * A(:, 1 : ccaDim) * codiag zeros(obj.dX, obj.D - ccaDim)];
    W2 = [yCov * B(:, 1 : ccaDim) * codiag zeros(obj.dY, obj.D - ccaDim)];
    obj.eX_W(:,:,k) = W1;
    obj.eY_W(:,:,k) = W2;

    obj.eX_psi(:,:,k) = 1*eye(obj.dX);
    obj.eY_psi(:,:,k) = 1*eye(obj.dY);

    %Dummy covariance matrices
    for j = 1 : obj.dX
        obj.covX_W(:,:,j,k) = eye(obj.D);
    end

    for j = 1 : obj.dY
        obj.covY_W(:,:,j,k) = eye(obj.D);
    end

    obj.covX_mu(:,:,k) = eye(obj.dX);
    obj.covY_mu(:,:,k) = eye(obj.dY);

    %Random latent coordinates:
    obj.cov_t(:,:,k) = eye(obj.D);
    obj.e_t(:,:,k) = randnorm(obj.N, zeros(obj.D,1),[], obj.cov_t(:,:,k));

    %ARD initialisation
    obj.e_alphaX = repmat(obj.ard_a / obj.ard_b, obj.D, obj.M);
    obj.e_alphaY = repmat(obj.ard_a / obj.ard_b, obj.D, obj.M);

end

obj.eX_gamma = repmat(obj.gammaX,obj.M,1);
obj.eY_gamma = repmat(obj.gammaY,obj.M,1);

%Latent scale initialisation
obj.e_u = ones(obj.N, obj.M);

end