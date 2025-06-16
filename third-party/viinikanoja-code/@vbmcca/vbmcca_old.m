% Variational mixture of robust CCAs
% Copyright (C) 2010 Jaakko Viinikanoja
%
% This library is free software; you can redistribute it and/or
% modify it under the terms of the GNU Lesser General Public
% License as published by the Free Software Foundation; either
% version 2.1 of the License, or (at your option) any later version.
%
% This library is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% Lesser General Public License for more details.
%
% You should have received a copy of the GNU Lesser General Public
% License along with this library; if not, write to the Free Software
% Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
%
%
% Depends on Lightspeed toolbox by Tom Minka
%

classdef vbmcca
    %VBMCCA Variational mixture of CCAs
    %   Documentation Check the source file & examples & function
    %   helps
    
    properties
        % General properties
        %--------------------------------------------------------------
        M;  %Cluster count
        dX; %Dimension of the X datavectors
        dY; %Dimension of the Y datavectors
        D;  %Dimension of the latent variables
        
        % Hyperparameters;
        betaX; % Mean parameter precision( inverse covariance scaling)
        betaY;
        
        ard_a; % ARD hyperparameters
        ard_b;
        
        gammaX % Data specific variation degrees of freedom (Wishart)
        gammaY;
        
        phiX;   % Data specific variation scale matrix
        phiY;
        
        mX; % Mean prior vectors in format dX X M
        mY;
        VX; % Projection matrix row prior vectors in format dX X D x M
        VY;
        
        % Additional parameters
        
        pi; % Cluster assignment probabalities
        nu; % t-distribution degrees of freedom
        
        % Learning associated properties
        %--------------------------------------------------------------
        
        N; %Data set size
        
        %Latent distribution statistics
        %------------------------------
        %Latent variable posterior statistics
        e_t; % D x N x M; D x N part describes latent variable means
        cov_t; %  D x D x M; Latent covariance without u dependency

        %Latent scales
        u_alpha;  %N x M; Gamma rv first argument
        u_beta; %N x M; Gamma rv second argument
        e_u;  %N x M; Scale first moments
        
        r;  %N x M; cluster assignments
       
        %Parameter posterior statistics
        %------------------------------
        eX_gamma;% M x 1; Wishart degrees of freedom
        eY_gamma;
        
        eX_psi; %dX x dX x M; Expected scales
        eY_psi;
        

        eX_mu; %dX x M; Mean posterior first moment
        eY_mu;
        
        covX_mu; % dX x dX x M; mean posterior covariance
        covY_mu;

        eX_W;   % d1 x D x M; Projection Row Means
        eY_W;   
        
        e_alphaX; % D X M; X-Projection matrix row precision
        alphaX_a;
        alphaX_b;
        
        e_alphaY; % D X M; Y-Projection matrix row precision
        alphaY_a;
        alphaY_b;

        covX_W; % D x D x dX x M; Covariance matrices
        covY_W;
        
        % Model Options
        %----------------------
        normalLatents; %Default is 0 for t-distribution, Switch to 1 to get Gaussian latents
        
        utDependency; % Assumed variational factorisation: 0 corresponds to u & t independent
                      % whilst 1 assumes u & t are dependent
    end
    
    methods
        function obj = vbmcca(M,dX,dY,D)
        % VBMCCA - Init the model and set the default hyperparameter
        % values.
        %
        % vbmcca(M,dX,dY,D) Arguments:
        % M - cluster count
        % dX - Data vector X dimension
        % dY - Data vector Y dimension
        % D - Latent variable dimension. 1 <= D <= min(dX,dY)
        %
        % You should set appropiates values for the hyperparameters after
        % initialising a vbmcca object.
        
            obj.M = M;
            obj.dX = dX;
            obj.dY = dY;
            obj.D = D;
            
            % Dummy hyperparameter initialisations:
            obj.betaX = 1;
            obj.betaY = 1;
        
            obj.ard_a = 0.1; 
            obj.ard_b = 0.1;
            
            obj.gammaX = dX + 1;
            obj.gammaY = dY + 1;
        
            obj.phiX = 1/100*eye(dX); %Inverse scale precision
            obj.phiY = 1/100*eye(dY);
        
            obj.mX = zeros(dX,M);
            obj.mY = zeros(dY,M);
       
            obj.VX = zeros(dX, D, M);
            obj.VY = zeros(dY, D, M);
       
            obj.pi = 1/M * ones(M,1);
            obj.nu = 10 * ones(M,1);
            
            obj.normalLatents = 0;
            obj.utDependency = 0;
        end
        
        function obj = initWithKMeans(obj, X, Y, iterations, retries)
        %INITWITHKMEANS Initialises the model parameters with k-means
        %clustering and CCA ML point estimates.
        %
        % obj = initWithKMeans(obj, X, Y, iterations, retries) 
        %
        % arguments:
        % X - dX X N datamatrix
        % Y - dY X N datamatrix
        % iterations - maximum number of k-means iterations, optional
        % retries - number of k-means repetitions, optional
        %
        % output:
        % obj - Model with initialised parameters.
       
            % Default parameter settings
            
            if(nargin <= 3)
                iterations = 150; 
            end
            
            if(nargin <= 4)
                retries = 10;
            end           
            
            % K-means clustering with retries
            
            data = [X' Y'];
            obj.N = min(cols(X), cols(Y));
            
            objectiveMin = Inf;
            indices = zeros(obj.N, 1);
            centers = zeros(obj.M, obj.dX + obj.dY);
            
            kMeansOptions =  statset('MaxIter',iterations,'Display','final');
         
            for i = 1 : retries
                [newIndices, newCenters,distance] = kmeans(data, obj.M, 'emptyaction', 'singleton','options', kMeansOptions);                
                    
                if( distance < objectiveMin)
                    centers = newCenters;
                    indices = newIndices;
                    objectiveMin = distance;
                end
                    
            end
            
            %Initialise posterior means from centers:
            obj.eX_mu = centers(:,1:obj.dX)';
            obj.eY_mu = centers(:,(obj.dX + 1):(obj.dX + obj.dY))';
            
            obj.r = zeros(obj.N, obj.M);
            
            %Initialise projections & data specific variance from prob cca
            %MLE results
            for k=1 : obj.M
                
                %clusterDataX = data(indices == k, 1:obj.dX);
                %clusterDataY = data(indices == k,(obj.dX + 1):(obj.dX + obj.dY));
                
                %Hard clustering
                obj.r(indices == k, k) = 1;
                
                %Cluster sizes
                obj.pi(k) = sum(indices == k) / obj.N;
                
                % Option for initialising with the classical CCA:
                % -------------------------------------------------
                %Data covariances
                %xCov = cov(clusterDataX);
                %yCov = cov(clusterDataY);
                
                %Normal CCA, Check Prob CCA paper for formula descriptions
                %[A,B,co] =  canoncorr(clusterDataX, clusterDataY);
                %ccaDim = obj.D;%min(rank(A'), rank(B'));
                %codiag = diag(sqrt( co(1 : ccaDim)));
                
                %obj.eX_W(:,:,k) = [xCov * A(:, 1 : ccaDim) * codiag zeros(obj.dX, obj.D - ccaDim)];
                %obj.eY_W(:,:,k) = [yCov * B(:, 1 : ccaDim) * codiag zeros(obj.dY, obj.D - ccaDim)];
                
                % Random initialisation:
                %-------------------------------------------------
                obj.eX_W(:,:,k) = rand(obj.dX, obj.D);
                obj.eY_W(:,:,k) = rand(obj.dY, obj.D);
                                
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
        
        function [obj,energy] = updateLatentPosterior(obj,X,Y, nuUpdates, initPhase)
        %UPDATELATENTPOSTERIOR - Updates the hidden variable posterior
        %statistics. For internal use.
        
            %Whether Hyperparameter nu is updated, default behaviour is
            %that updating is enabled
            if nargin < 3
                nuUpdates = 1;
            end
        
            % Auxiliary variables for evaluation
            %------------------------------------------------------------
            centeredX = zeros(obj.dX,obj.N,obj.M);
            centeredY = zeros(obj.dY,obj.N,obj.M);
            
            eX_WpsiW = zeros(obj.D, obj.D, obj.M);
            eY_WpsiW = zeros(obj.D, obj.D, obj.M);
            
            eX_psiDiagCovW = zeros(obj.D,obj.D,obj.M);
            eY_psiDiagCovW = zeros(obj.D,obj.D,obj.M);
            
            e_logu = zeros(obj.N, obj.M);
            
            DKL_t = zeros(obj.N, obj.M);
            DKL_u = zeros(obj.N, obj.M);
            %DKL_p = zeros(obj.N, obj.M);
            
            % Precentered data matrices
            for k = 1 : obj.M
                centeredX(:,:,k) = X - repmat(obj.eX_mu(:,k),1,obj.N);
                centeredY(:,:,k) = Y - repmat(obj.eY_mu(:,k),1,obj.N);
            end  
            
            % Precalculation of auxiliary parameter matrices
            for k = 1 : obj.M
                
                % Assembly auxiliary matrices
                for j = 1 : obj.dX
                    eX_psiDiagCovW(:,:,k) = eX_psiDiagCovW(:,:,k) + obj.eX_psi(j,j,k) * obj.covX_W(:,:,j,k);
                end
                
                for j = 1 : obj.dY
                    eY_psiDiagCovW(:,:,k) = eY_psiDiagCovW(:,:,k) + obj.eY_psi(j,j,k) * obj.covY_W(:,:,j,k);
                end
                
                eX_WpsiW(:,:,k) = obj.eX_W(:,:,k)' * obj.eX_psi(:,:,k) * obj.eX_W(:,:,k) + eX_psiDiagCovW(:,:,k);
                eY_WpsiW(:,:,k) = obj.eY_W(:,:,k)' * obj.eY_psi(:,:,k) * obj.eY_W(:,:,k) + eY_psiDiagCovW(:,:,k);
                            
            end
            
            % t-latent Posterior
            %-------------------------------------------------------------
      
            %
            % t-means & without u covariances \Sigma_{t_k}
            %
            
            for k = 1 : obj.M    
                
                obj.cov_t(:,:,k) = inv(eye(obj.D) + eX_WpsiW(:,:,k) + eY_WpsiW(:,:,k));
                
                obj.e_t(:,:,k) = obj.cov_t(:,:,k) * ( obj.eX_W(:,:,k)' * obj.eX_psi(:,:,k) * ...
                    centeredX(:,:,k) + obj.eY_W(:,:,k)' * obj.eY_psi(:,:,k) * centeredY(:,:,k));
            end
            
            % Gamma Posterior
            %-------------------------------------------------------------
            
            if obj.normalLatents == 0
        
                %Variational approximation is of form q(t|u,z) & q(u|z)
                if obj.utDependency == 1
                    %Posterior alphas
                    obj.u_alpha = repmat(0.5*(obj.nu' + (obj.dX + obj.dY )), obj.N,1);
            
                    % Posterior betas
                    for k = 1 : obj.M
                    
                        xContrib = 0.5 * col_sum((centeredX(:,:,k) - obj.eX_W(:,:,k)*obj.e_t(:,:,k)) .* (obj.eX_psi(:,:,k)* ...
                            (centeredX(:,:,k) - obj.eX_W(:,:,k)*obj.e_t(:,:,k))))' ...
                            + 0.5 * trace(obj.eX_psi(:,:,k)*obj.covX_mu(:,:,k)) * ones(obj.N,1) ...
                            + 0.5 * col_sum(obj.e_t(:,:,k) .* (eX_psiDiagCovW(:,:,k) *obj.e_t(:,:,k)))';
                                        
                        yContrib = 0.5 * col_sum((centeredY(:,:,k) - obj.eY_W(:,:,k)*obj.e_t(:,:,k)) .* (obj.eY_psi(:,:,k)* ...
                            (centeredY(:,:,k) - obj.eY_W(:,:,k)*obj.e_t(:,:,k))))' ...
                            + 0.5 * trace(obj.eY_psi(:,:,k)*obj.covY_mu(:,:,k)) * ones(obj.N,1) ...
                            + 0.5 * col_sum(obj.e_t(:,:,k) .* (eY_psiDiagCovW(:,:,k)*obj.e_t(:,:,k)))';
                    
                        obj.u_beta(:, k) = xContrib + yContrib + 0.5*col_sum((obj.e_t(:,:,k).*obj.e_t(:,:,k)))' ...
                            + 0.5*obj.nu(k)*ones(obj.N,1);
                      
                    end
                    
                %Variational approximation is of form q(t|z) & q(u|z)
                else
                    obj.u_alpha = repmat(0.5*(obj.nu' + (obj.dX + obj.dY + obj.D)), obj.N,1);
            
                    % Posterior betas
                    for k = 1 : obj.M
                    
                        xContrib = 0.5 * col_sum((centeredX(:,:,k) - obj.eX_W(:,:,k)*obj.e_t(:,:,k)) .* (obj.eX_psi(:,:,k)* ...
                            (centeredX(:,:,k) - obj.eX_W(:,:,k)*obj.e_t(:,:,k))))' ...
                            + 0.5 * trace(obj.eX_psi(:,:,k)*obj.covX_mu(:,:,k)) * ones(obj.N,1) ...
                            + 0.5 * col_sum(obj.e_t(:,:,k) .* (eX_psiDiagCovW(:,:,k) *obj.e_t(:,:,k)))' ...
                            + 0.5 * trace(obj.cov_t(:,:,k) * eX_WpsiW(:,:,k)) * (ones(obj.N,1) ./obj.e_u(:,k));
                                        
                        yContrib = 0.5 * col_sum((centeredY(:,:,k) - obj.eY_W(:,:,k)*obj.e_t(:,:,k)) .* (obj.eY_psi(:,:,k)* ...
                            (centeredY(:,:,k) - obj.eY_W(:,:,k)*obj.e_t(:,:,k))))' ...
                            + 0.5 * trace(obj.eY_psi(:,:,k)*obj.covY_mu(:,:,k)) * ones(obj.N,1) ...
                            + 0.5 * col_sum(obj.e_t(:,:,k) .* (eY_psiDiagCovW(:,:,k)*obj.e_t(:,:,k)))' ...
                            + 0.5 * trace(obj.cov_t(:,:,k) * eY_WpsiW(:,:,k)) * (ones(obj.N,1) ./obj.e_u(:,k));
                    
                        obj.u_beta(:, k) = xContrib + yContrib + 0.5*col_sum((obj.e_t(:,:,k).*obj.e_t(:,:,k)))' ...
                            + 0.5 * trace(obj.cov_t(:,:,k)) * (ones(obj.N,1) ./obj.e_u(:,k)) + 0.5*obj.nu(k)*ones(obj.N,1);
                      
                    end
                end
            
                %First moment of gamma distribution
                obj.e_u = obj.u_alpha ./ obj.u_beta;
                %Expectation of the logarithm
                e_logu = digamma(obj.u_alpha) - log(obj.u_beta);
            else
            
                %Force normality
                obj.e_u = ones(obj.N,obj.M);
                %logarithms are initialised to zero by default
                
            end
            
            % Estimate degrees of freedom by calculating the point estimate
            % of nu. 
            %------------------------------------------------------------
            if obj.normalLatents == 0
            
                degreesOfFreedom = zeros(obj.M, 1);
                %Numerical optimization for nu:s
                for k = 1 : obj.M
                    %Pre-evalution for a function
                
                    constant = 1/(obj.N*obj.pi(k)) *(obj.r(:,k)' * (e_logu(:,k) - obj.e_u(:,k) ));
                
                    %Goal function
                    goalFunction = @(x) (1 + log(x/2) - digamma(x/2) + constant);
                
                    %Start iteration from previous value
                    degreesOfFreedom(k) = fsolve(goalFunction, obj.nu(k), optimset('Display','off'));
                end
                
                %Point estimate update
                if nuUpdates == 1
                    obj.nu = degreesOfFreedom;
                end
                
                disp(strcat('Estimated degrees of freedom: ',num2str(degreesOfFreedom')))
            end
                                
            %
            % Cluster responsibilities
            %-------------------------------------------------------------
            
            % KL-divergences needed for cluster responsibilities
            
            for k = 1: obj.M
            
                %This includes evaluation of expectation wrt. q(u_nk | z_nk) 
                DKL_t(:,k) = 0.5*(repmat(trace(obj.cov_t(:,:,k)) - logdet(obj.cov_t(:,:,k)) - obj.D ,obj.N,1) + ...
                    (col_sum(obj.e_t(:,:,k) .* obj.e_t(:,:,k)))' .* obj.e_u(:,k));
                
                %Alpha & beta parameter for p-gamma distribution
                p_ab = repmat(0.5 * obj.nu',obj.N,1);
                
                if obj.normalLatents == 0
                
                    DKL_u = gammaln(p_ab) - gammaln(obj.u_alpha) + obj.u_alpha .* log(obj.u_beta) - ...
                        p_ab .* log(p_ab) + (obj.u_alpha - p_ab).*(digamma(obj.u_alpha) - log(obj.u_beta)) + ...
                        obj.u_alpha .* ((p_ab - obj.u_beta) ./ obj.u_beta);
                end
                %"else" Otherwise keep zero entries from initialisation 
            end
            
            elnX_p = zeros(obj.N, obj.M);
            elnY_p = zeros(obj.N, obj.M);
            
            for k = 1 : obj.M
                
                % Wishart det logarithm expectations
                elnX_detPsi = vbmcca.eLogDetPsi(obj.eX_psi(:,:,k), obj.eX_gamma(k),obj.dX);
                elnY_detPsi = vbmcca.eLogDetPsi(obj.eY_psi(:,:,k), obj.eY_gamma(k),obj.dY);
               
                % likelihoods       
                
                elnX_p(:,k) = -0.5 * obj.dX * log(2*pi) * ones(obj.N,1) ...
                    +0.5*(repmat(elnX_detPsi,obj.N,1) + obj.dX * e_logu(:,k)) ...
                    -0.5 * col_sum((centeredX(:,:,k) - obj.eX_W(:,:,k)*obj.e_t(:,:,k)) .* (obj.eX_psi(:,:,k)* ...
                    (centeredX(:,:,k) - obj.eX_W(:,:,k)*obj.e_t(:,:,k))))' .* obj.e_u(:,k) ...
                    -0.5 * trace(obj.eX_psi(:,:,k)*obj.covX_mu(:,:,k)) * ones(obj.N,1) .* obj.e_u(:,k)  ...
                    -0.5 * trace(obj.cov_t(:,:,k) * eX_WpsiW(:,:,k)) * ones(obj.N,1) ...
                    -0.5 * col_sum(obj.e_t(:,:,k) .* (eX_psiDiagCovW(:,:,k)*obj.e_t(:,:,k)))' .* obj.e_u(:,k);
               
                
                elnY_p(:,k) = -0.5 * obj.dY * log(2*pi) * ones(obj.N,1) ...
                    +0.5*(repmat(elnY_detPsi,obj.N,1) + obj.dY* e_logu(:,k)) ...
                    -0.5 * col_sum((centeredY(:,:,k) - obj.eY_W(:,:,k)*obj.e_t(:,:,k)) .* (obj.eY_psi(:,:,k)* ...
                    (centeredY(:,:,k) - obj.eY_W(:,:,k)*obj.e_t(:,:,k))))' .* obj.e_u(:,k) ...
                    -0.5 * trace(obj.eY_psi(:,:,k)*obj.covY_mu(:,:,k)) * ones(obj.N,1) .* obj.e_u(:,k)  ...
                    -0.5 * trace(obj.cov_t(:,:,k) * eY_WpsiW(:,:,k)) * ones(obj.N,1) ...
                    -0.5 * col_sum(obj.e_t(:,:,k) .* (eY_psiDiagCovW(:,:,k)*obj.e_t(:,:,k)))' .* obj.e_u(:,k);
            
            end
            
            % Calculate rhos / rs
            rhos = exp( repmat(log(obj.pi'),obj.N,1) - DKL_u - DKL_t + elnX_p + elnY_p);
            obj.r = scale_rows(rhos, 1 ./ row_sum(rhos));
            
            % Play around with possibly very bad starting values
            if (initPhase == 1)
                obj.r(isnan(obj.r)) = 1/obj.M;
                obj.r(isinf(obj.r)) = 1/obj.M;
            end
           
            % Free energy (latent variable part + likelihoods)
            %-------------------------------------------------------------
            
            %p KL-divergence for p
            DKL_p = obj.r .* log((obj.r)./ repmat(obj.pi',obj.N,1));
            %Fix zero limits
            indices = isnan(DKL_p);
            DKL_p(indices) = 0;
            indices = isinf(DKL_p);
            DKL_p(indices) = 0;
            
            energy = sum(row_sum(obj.r .* (elnX_p + elnY_p - DKL_t - DKL_u))) - sum(row_sum(DKL_p));
            
            % Type 2 ML estimates for pi:s
            %--------------------------------------------------------------
            
            %pi
            obj.pi = 1/obj.N * col_sum(obj.r)';
            
            
        end
        
        function statistics = inferLatentsX(obj, X, maxuIter, uTolerance)
            %INFERLATENTSX - Infer the posterior distributions for new
            %samples using the given model(obj), each column of the data
            %matrix X should correspond to a datavector in view X. 
            %
            % statistics = inferLatentsX(obj, X, maxuIter, uTolerance)
            % Parameters:
            % X - the dataset for which latents are inferred. Dim(X) = dX X newN
            % maxuIter - The maximum number of iterations in u optimisations,
            % optional
            % uTolerance - Alternative termination condition for u optimization
            % iteration termination, optional
             
            % Default parameter values:  ( maxIter = 100, uTolerance =
            % 10^(-4) )
            % 
            % output:
            % statistics - Inferred latenStatistics structure
            
            %-----------------------------------------------------------
            if(nargin == 2)
                maxuIter = 100;
                uTolerance = 10^(-4);
            end
          
            
            % Auxiliary variables for evaluation
            %------------------------------------------------------------
            
            newdataN = cols(X);
            
            centeredX = zeros(obj.dX, newdataN, obj.M);
            eX_WpsiW = zeros(obj.D, obj.D, obj.M);
            eX_psiDiagCovW = zeros(obj.D,obj.D,obj.M);
  
            e_newdata_u = ones(newdataN, obj.M);
            e_newdata_logu = zeros(newdataN, obj.M);
            u_newdata_alpha = zeros(newdataN, obj.M);
            u_newdata_beta = zeros(newdataN, obj.M);
            
            e_newdata_t = zeros(obj.D, newdataN, obj.M);
            cov_newdata_t = zeros(obj.D, obj.D, obj.M);
            
            DKL_t = zeros(newdataN, obj.M);
            DKL_u = zeros(newdataN, obj.M);
            %DKL_p = zeros(newdataN, obj.M);
            
            % Precentered data matrices
            for k = 1 : obj.M
                centeredX(:,:,k) = X - repmat(obj.eX_mu(:,k), 1, newdataN);
            end  
            
            % Precalculation of auxiliary parameter matrices
            for k = 1 : obj.M
                
                % Assembly auxiliary matrices
                for j = 1 : obj.dX
                    eX_psiDiagCovW(:,:,k) = eX_psiDiagCovW(:,:,k) + obj.eX_psi(j,j,k) * obj.covX_W(:,:,j,k);
                end
                eX_WpsiW(:,:,k) = obj.eX_W(:,:,k)' * obj.eX_psi(:,:,k) * obj.eX_W(:,:,k) + eX_psiDiagCovW(:,:,k);
                            
            end
            
            % t-latent Posterior
            %-------------------------------------------------------------
      
            %
            % t-means & without u covariances \Sigma_{t_k}
            %
            
            for k = 1 : obj.M    
                
                cov_newdata_t(:,:,k) = inv(eye(obj.D) + eX_WpsiW(:,:,k));
                
                e_newdata_t(:,:,k) = cov_newdata_t(:,:,k) * ( obj.eX_W(:,:,k)' * obj.eX_psi(:,:,k) * ...
                    centeredX(:,:,k));
            end
            
            % Gamma Posterior
            %-------------------------------------------------------------
            
            if obj.normalLatents == 0
        
                 if obj.utDependency == 0
                    %Posterior alphas
                    u_newdata_alpha = repmat(0.5*(obj.nu' + obj.dX + obj.D), newdataN,1);
                    previous_e_u = e_newdata_u;
                     
                    %Iterative optimization
                    for j = 1 : maxuIter
                        
                        % Posterior betas
                        for k = 1 : obj.M
                    
                            xContrib = 0.5 * col_sum((centeredX(:,:,k) - obj.eX_W(:,:,k)*e_newdata_t(:,:,k)) .* (obj.eX_psi(:,:,k)* ...
                                (centeredX(:,:,k) - obj.eX_W(:,:,k)*e_newdata_t(:,:,k))))' ...
                                + 0.5 * trace(obj.eX_psi(:,:,k)*obj.covX_mu(:,:,k)) * ones(newdataN,1) ...
                                + 0.5 * col_sum(e_newdata_t(:,:,k) .* (eX_psiDiagCovW(:,:,k) *e_newdata_t(:,:,k)))'...
                                + 0.5 * trace(cov_newdata_t(:,:,k) * eX_WpsiW(:,:,k)) * (ones(newdataN,1) ./ e_newdata_u(:,k));
                                          
                            u_newdata_beta(:, k) = xContrib + 0.5*col_sum((e_newdata_t(:,:,k).*e_newdata_t(:,:,k)))' ...
                                + 0.5 * trace(cov_newdata_t(:,:,k)) * (ones(newdataN,1) ./ e_newdata_u(:,k)) + 0.5*obj.nu(k)*ones(newdataN,1);
                        end
                      
                        %First moment of gamma distribution
                        e_newdata_u = u_newdata_alpha ./ u_newdata_beta;
                        %Expectation of the logarithm
                        e_newdata_logu = digamma(u_newdata_alpha) - log(u_newdata_beta);
                        
                        %Convergence check
                        if( (norm(previous_e_u - e_newdata_u) / rows(e_newdata_u)) < uTolerance)
                            %disp(strcat('U optimisation terminated after: ',num2str(j),' iterations'));
                            break;
                        end
                        previous_e_u = e_newdata_u;
                    end
                 else
                    %Gaussian-Gamma nonfactorised posterior
                     %Posterior alphas
                    u_newdata_alpha = repmat(0.5*(obj.nu' + obj.dX ), newdataN,1);
            
                    % Posterior betas
                    for k = 1 : obj.M
                    
                        xContrib = 0.5 * col_sum((centeredX(:,:,k) - obj.eX_W(:,:,k)*e_newdata_t(:,:,k)) .* (obj.eX_psi(:,:,k)* ...
                            (centeredX(:,:,k) - obj.eX_W(:,:,k)*e_newdata_t(:,:,k))))' ...
                            + 0.5 * trace(obj.eX_psi(:,:,k)*obj.covX_mu(:,:,k)) * ones(newdataN,1) ...
                            + 0.5 * col_sum(e_newdata_t(:,:,k) .* (eX_psiDiagCovW(:,:,k) *e_newdata_t(:,:,k)))';                                      
                    
                        u_newdata_beta(:, k) = xContrib + 0.5*col_sum((e_newdata_t(:,:,k).*e_newdata_t(:,:,k)))' ...
                            + 0.5*obj.nu(k)*ones(newdataN,1);
                      
                    end
                    e_newdata_u = u_newdata_alpha ./ u_newdata_beta;
                    e_newdata_logu = digamma(u_newdata_alpha) - log(u_newdata_beta);
                end          
            end
                              
            %
            % Cluster responsibilities
            %-------------------------------------------------------------
    
            for k = 1: obj.M
            
                %This includes evaluation of expectation wrt. q(u_nk | z_nk) 
                DKL_t(:,k) = 0.5*(repmat(trace(cov_newdata_t(:,:,k)) - logdet(cov_newdata_t(:,:,k)) - obj.D, newdataN, 1) + ...
                    (col_sum(e_newdata_t(:,:,k) .* e_newdata_t(:,:,k)))' .* e_newdata_u(:,k));
                
                %Alpha & beta parameter for p-gamma distribution
                p_ab = repmat(0.5 * obj.nu', newdataN,1);
                
                if obj.normalLatents == 0
                
                    DKL_u = gammaln(p_ab) - gammaln(u_newdata_alpha) + u_newdata_alpha .* log(u_newdata_beta) - ...
                        p_ab .* log(p_ab) + (u_newdata_alpha - p_ab).*(digamma(u_newdata_alpha) - log(u_newdata_beta)) + ...
                        u_newdata_alpha .* ((p_ab - u_newdata_beta) ./ u_newdata_beta);
                end
                %Otherwise keep zero entries from initialisation    
            end
            
            elnX_p = zeros(newdataN, obj.M);
            
            for k = 1 : obj.M
                
                % Wishart det logarithm expectations
                elnX_detPsi = vbmcca.eLogDetPsi(obj.eX_psi(:,:,k), obj.eX_gamma(k),obj.dX);
               
                % likelihoods       
                
                elnX_p(:,k) = -0.5 * obj.dX * log(2*pi) * ones(newdataN, 1) ...
                    +0.5*(repmat(elnX_detPsi, newdataN, 1) + obj.dX * e_newdata_logu(:,k)) ...
                    -0.5 * col_sum((centeredX(:,:,k) - obj.eX_W(:,:,k)* e_newdata_t(:,:,k)) .* (obj.eX_psi(:,:,k)* ...
                    (centeredX(:,:,k) - obj.eX_W(:,:,k)*e_newdata_t(:,:,k))))' .* e_newdata_u(:,k) ...
                    -0.5 * trace(obj.eX_psi(:,:,k)*obj.covX_mu(:,:,k)) * ones(newdataN,1) .* e_newdata_u(:,k)  ...
                    -0.5 * trace(cov_newdata_t(:,:,k) * eX_WpsiW(:,:,k)) * ones(newdataN, 1) ...
                    -0.5 * col_sum(e_newdata_t(:,:,k) .* (eX_psiDiagCovW(:,:,k)*e_newdata_t(:,:,k)))' .* e_newdata_u(:,k);
              
            end
            
            % Calculate rhos / rs
            rhos = exp( repmat(log(obj.pi'), newdataN, 1) - DKL_u - DKL_t + elnX_p);
            newdata_r = scale_rows(rhos, 1 ./ row_sum(rhos));
            
            % Save comparison statistics
            statistics = obj.initLatentStatistics(obj.D, newdataN, obj.M);
            statistics.e_t = e_newdata_t;
            statistics.cov_t = cov_newdata_t;
            statistics.u_alpha = u_newdata_alpha;
            statistics.u_beta = u_newdata_beta;
            statistics.cov_t_weights = 1 ./ e_newdata_u;
            statistics.r = newdata_r;       
        end
        
        
       function statistics = inferLatentsY(obj, Y, maxuIter, uTolerance)
            %INFERLATENTSY - Infer the posterior distributions for new
            %samples using the given model(obj), each column of the data
            %matrix Y should correspond to a datavector in view Y. 
            %
            % statistics = inferLatentsX(obj, Y, maxuIter, uTolerance)
            % Parameters:
            % Y - the dataset for which latents are inferred. Dim(Y) = dY X newN
            % maxuIter - The maximum number of iterations in u optimisations,
            % optional
            % uTolerance - Alternative termination condition for u optimization
            % iteration termination, optional
             
            % Default parameter values:  ( maxIter = 100, uTolerance =
            % 10^(-4) )
            % 
            % output:
            % statistics - Inferred latenStatistics structure
             
            %-----------------------------------------------------------
            if(nargin == 2)
                maxuIter = 100;
                uTolerance = 10^(-4);
            end
          
            
            % Auxiliary variables for evaluation
            %------------------------------------------------------------
            
            newdataN = cols(Y);
            
            centeredY = zeros(obj.dX, newdataN, obj.M);
            eY_WpsiW = zeros(obj.D, obj.D, obj.M);
            eY_psiDiagCovW = zeros(obj.D,obj.D,obj.M);
  
            e_newdata_u = ones(newdataN, obj.M);
            e_newdata_logu = zeros(newdataN, obj.M);
            u_newdata_alpha = zeros(newdataN, obj.M);
            u_newdata_beta = zeros(newdataN, obj.M);
            
            e_newdata_t = zeros(obj.D, newdataN, obj.M);
            cov_newdata_t = zeros(obj.D, obj.D, obj.M);
            
            DKL_t = zeros(newdataN, obj.M);
            DKL_u = zeros(newdataN, obj.M);
            %DKL_p = zeros(newdataN, obj.M);
            
            % Precentered data matrices
            for k = 1 : obj.M
                centeredY(:,:,k) = Y - repmat(obj.eY_mu(:,k), 1, newdataN);
            end  
            
            % Precalculation of auxiliary parameter matrices
            for k = 1 : obj.M
                
                % Assembly auxiliary matrices
                for j = 1 : obj.dX
                    eY_psiDiagCovW(:,:,k) = eY_psiDiagCovW(:,:,k) + obj.eY_psi(j,j,k) * obj.covY_W(:,:,j,k);
                end
                eY_WpsiW(:,:,k) = obj.eY_W(:,:,k)' * obj.eY_psi(:,:,k) * obj.eY_W(:,:,k) + eY_psiDiagCovW(:,:,k);
                            
            end
            
            % t-latent Posterior
            %-------------------------------------------------------------
      
            % t-means & without u covariances \Sigma_{t_k}
            for k = 1 : obj.M    
                
                cov_newdata_t(:,:,k) = inv(eye(obj.D) + eY_WpsiW(:,:,k));
                
                e_newdata_t(:,:,k) = cov_newdata_t(:,:,k) * ( obj.eY_W(:,:,k)' * obj.eY_psi(:,:,k) * ...
                    centeredY(:,:,k));
            end
            
            % Gamma Posterior
            %-------------------------------------------------------------
            
            if obj.normalLatents == 0
                     
                if obj.utDependency == 0
                    %Posterior alphas
                    u_newdata_alpha = repmat(0.5*(obj.nu' + obj.dY + obj.D), newdataN,1);
                    previous_e_u = e_newdata_u;
                     
                    %Iterative optimization
                    for j = 1 : maxuIter
                        
                        % Posterior betas
                        for k = 1 : obj.M
                    
                            yContrib = 0.5 * col_sum((centeredY(:,:,k) - obj.eY_W(:,:,k)*e_newdata_t(:,:,k)) .* (obj.eY_psi(:,:,k)* ...
                                (centeredY(:,:,k) - obj.eY_W(:,:,k)*e_newdata_t(:,:,k))))' ...
                                + 0.5 * trace(obj.eY_psi(:,:,k)*obj.covY_mu(:,:,k)) * ones(newdataN,1) ...
                                + 0.5 * col_sum(e_newdata_t(:,:,k) .* (eY_psiDiagCovW(:,:,k) *e_newdata_t(:,:,k)))'...
                                + 0.5 * trace(cov_newdata_t(:,:,k) * eY_WpsiW(:,:,k)) * (ones(newdataN,1) ./ e_newdata_u(:,k));
                                          
                            u_newdata_beta(:, k) = yContrib + 0.5*col_sum((e_newdata_t(:,:,k).*e_newdata_t(:,:,k)))' ...
                            + 0.5 * trace(cov_newdata_t(:,:,k)) * (ones(newdataN,1) ./ e_newdata_u(:,k)) + 0.5*obj.nu(k)*ones(newdataN,1);
                        end
                      
                        %First moment of gamma distribution
                        e_newdata_u = u_newdata_alpha ./ u_newdata_beta;
                        %Expectation of the logarithm
                        e_newdata_logu = digamma(u_newdata_alpha) - log(u_newdata_beta);
                        
                        %Convergence check
                        if( (norm(previous_e_u - e_newdata_u) / rows(e_newdata_u)) < uTolerance)
                            %disp(strcat('U optimisation terminated after: ',num2str(j),' iterations'));
                            break;
                        end
                        previous_e_u = e_newdata_u;
                    end
                else
                    %Posterior alphas
                    u_newdata_alpha = repmat(0.5*(obj.nu'+ obj.dY), newdataN,1);
            
                    % Posterior betas
                    for k = 1 : obj.M
                    
            
                                        
                        yContrib = 0.5 * col_sum((centeredY(:,:,k) - obj.eY_W(:,:,k)*e_newdata_t(:,:,k)) .* (obj.eY_psi(:,:,k)* ...
                            (centeredY(:,:,k) - obj.eY_W(:,:,k)*e_newdata_t(:,:,k))))' ...
                            + 0.5 * trace(obj.eY_psi(:,:,k)*obj.covY_mu(:,:,k)) * ones(newdataN,1) ...
                            + 0.5 * col_sum(e_newdata_t(:,:,k) .* (eY_psiDiagCovW(:,:,k)*e_newdata_t(:,:,k)))';
                    
                        u_newdata_beta(:, k) = yContrib + 0.5*col_sum((e_newdata_t(:,:,k).*e_newdata_t(:,:,k)))' ...
                            + 0.5*obj.nu(k)*ones(newdataN,1);
                      
                    end
                    
                    %First moment of gamma distribution
                    e_newdata_u = u_newdata_alpha ./ u_newdata_beta;
                    %Expectation of the logarithm
                    e_newdata_logu = digamma(u_newdata_alpha) - log(u_newdata_beta);
                end
            end
                              
            %
            % Cluster responsibilities
            %-------------------------------------------------------------
    
            for k = 1: obj.M
            
                %This includes evaluation of expectation wrt. q(u_nk | z_nk) 
                DKL_t(:,k) = 0.5*(repmat(trace(cov_newdata_t(:,:,k)) - logdet(cov_newdata_t(:,:,k)) - obj.D, newdataN, 1) + ...
                    (col_sum(e_newdata_t(:,:,k) .* e_newdata_t(:,:,k)))' .* e_newdata_u(:,k));
                
                %Alpha & beta parameter for p-gamma distribution
                p_ab = repmat(0.5 * obj.nu', newdataN,1);
                
                if obj.normalLatents == 0
                
                    DKL_u = gammaln(p_ab) - gammaln(u_newdata_alpha) + u_newdata_alpha .* log(u_newdata_beta) - ...
                        p_ab .* log(p_ab) + (u_newdata_alpha - p_ab).*(digamma(u_newdata_alpha) - log(u_newdata_beta)) + ...
                        u_newdata_alpha .* ((p_ab - u_newdata_beta) ./ u_newdata_beta);
                end
                %Otherwise keep zero entries from initialisation    
            end
            
            elnY_p = zeros(newdataN, obj.M);
            
            for k = 1 : obj.M
                
                % Wishart det logarithm expectations
                elnY_detPsi = vbmcca.eLogDetPsi(obj.eX_psi(:,:,k), obj.eX_gamma(k),obj.dX);
               
                % likelihoods       
                
                elnY_p(:,k) = -0.5 * obj.dY * log(2*pi) * ones(newdataN, 1) ...
                    +0.5*(repmat(elnY_detPsi, newdataN, 1) + obj.dY * e_newdata_logu(:,k)) ...
                    -0.5 * col_sum((centeredY(:,:,k) - obj.eY_W(:,:,k)* e_newdata_t(:,:,k)) .* (obj.eY_psi(:,:,k)* ...
                    (centeredY(:,:,k) - obj.eY_W(:,:,k)*e_newdata_t(:,:,k))))' .* e_newdata_u(:,k) ...
                    -0.5 * trace(obj.eY_psi(:,:,k)*obj.covY_mu(:,:,k)) * ones(newdataN,1) .* e_newdata_u(:,k)  ...
                    -0.5 * trace(cov_newdata_t(:,:,k) * eY_WpsiW(:,:,k)) * ones(newdataN, 1) ...
                    -0.5 * col_sum(e_newdata_t(:,:,k) .* (eY_psiDiagCovW(:,:,k)*e_newdata_t(:,:,k)))' .* e_newdata_u(:,k);
              
            end
            
            % Calculate rhos / rs
            rhos = exp( repmat(log(obj.pi'), newdataN, 1) - DKL_u - DKL_t + elnY_p);
            newdata_r = scale_rows(rhos, 1 ./ row_sum(rhos));
            
            % Save the sufficient latent statistics
            statistics = obj.initLatentStatistics(obj.D, newdataN, obj.M);
            statistics.e_t = e_newdata_t;
            statistics.cov_t = cov_newdata_t;
            statistics.u_alpha = u_newdata_alpha;
            statistics.u_beta = u_newdata_beta;
            statistics.cov_t_weights = 1 ./ e_newdata_u;
            statistics.r = newdata_r;       
       end
        
       function [eY_mu_pred,eY_cov_pred,latentStatistics] = inferConditionalYMean(obj, X)
       % InferConditionalYMean - Calculates the conditional expectation of
       % Y when X is given
       %
       % [eY_mu_pred,eY_cov_pred,latentStatistics] =
       % inferConditionalYMean(obj, X)
       %
       % Arguments:
       % obj - X is a dX × newN data-matrix for new data in view X.
       %
       % Output:
       % eY_mu_pred - Conditional expectation wrt. posterior E[E[Y|X]]
       % eY_cov_pred - Conditional covariance wrt. posterior 
       % latentStatistics - Inferred latent variable statistics
       
            latentStatistics = obj.inferLatentsX(X);
            eY_mu_pred = zeros(obj.dY, cols(X));
            eY_cov_pred = zeros(obj.dY,obj.dY, cols(X));
            
            %Means
            for k = 1 : obj.M
                eY_mu_pred = eY_mu_pred + scale_cols( obj.eY_W(:,:,k) * latentStatistics.e_t(:,:,k) ...
                    + repmat(obj.eY_mu(:,k),1,cols(X)), latentStatistics.r(:,k));
            end
           
            %Inverse Wishart expectation
            e_inv_psi = zeros(obj.dY, obj.dY, obj.M);
            %Inverse gamma expectation
            if obj.normalLatents == 0
                e_inv_u = 1 ./ (latentStatistics.u_beta .* (latentStatistics.u_alpha -1));
            else
                e_inv_u = ones(cols(X),obj.M);
            end
                
            for k = 1:obj.M
                e_inv_psi(:,:,k) = inv(obj.eY_psi(:,:,k) / obj.eY_gamma(k)) / (obj.eY_gamma(k) - obj.dY - 1);
            end
            
            
            for n = 1: cols(X)
                          
                for k = 1: obj.M
                    
                    %Auxiliary matrices
                    e_tt = latentStatistics.e_t(:,n,k)*latentStatistics.e_t(:,n,k)' + latentStatistics.cov_t(:,:,k) * latentStatistics.cov_t_weights(n,k);
                    diagTraces = zeros(obj.dY,obj.dY);
                    
                    for j = 1 : obj.dY
                        diagTraces(j,j) = trace(e_tt * obj.covY_W(:,:,j,k));
                    end
                    
                    
                    %cluster specific contributions
                    eY_cov_pred(:, :, n) = eY_cov_pred(:, :, n) + latentStatistics.r(n,k) * (e_inv_u(n,k) * e_inv_psi(:,:,k) ...
                        + obj.eY_W(:,:,k) * latentStatistics.e_t(:,n,k) * obj.eY_mu(:,k)' ...
                        + obj.eY_mu(:,k) * latentStatistics.e_t(:,n,k)' * obj.eY_W(:,:,k)' ...
                        + obj.covY_mu(:,:,k) + obj.eY_mu(:,k)*obj.eY_mu(:,k)' ...
                        + obj.eY_W(:,:,k) *(e_tt) * obj.eY_W(:,:,k)' ...
                        + diagTraces ); 
                end
                %Covariances
                eY_cov_pred(:, :, n) = eY_cov_pred(:, :, n) - eY_mu_pred(:,n) * eY_mu_pred(:,n)';
            end
            
       end
       
       function [eX_mu_pred,eX_cov_pred] = inferConditionalXMean(obj, Y)
       % InferConditionalXMean - Calculates the conditional expectation of
       % X when Y is given
       %
       % [eX_mu_pred,eX_cov_pred] = inferConditionalXMean(obj, Y)
       %
       % Arguments:
       % obj - Y is a dY × newN data-matrix for new data in view Y.
       %
       % Output:
       % eX_mu_pred - Conditional expectation wrt. posterior E[E[X|Y]]
       % eX_cov_pred - Conditional covariance wrt. posterior 
       % latentStatistics - Inferred latent variable statistics
       
            latentStatistics = obj.inferLatentsY(Y);
            eX_mu_pred = zeros(obj.dX, cols(Y));
            eX_cov_pred = zeros(obj.dX,obj.dX, cols(Y));

            
            %Sweep over the clusters
            for k = 1 : obj.M
                eX_mu_pred = eX_mu_pred + scale_cols( obj.eX_W(:,:,k) * latentStatistics.e_t(:,:,k) ...
                    + repmat(obj.eX_mu(:,k),1,cols(Y)), latentStatistics.r(:,k));
            end
            
            %Inverse Wishart expectation
            e_inv_psi = zeros(obj.dX, obj.dX, obj.M);
            %Inverse gamma expectation
            if obj.normalLatents == 0
                 e_inv_u = 1 ./ (latentStatistics.u_beta .* (latentStatistics.u_alpha -1));
            else
                e_inv_u = ones(cols(Y),obj.M);
            end
                
            for k = 1:obj.M
                e_inv_psi(:,:,k) = inv(obj.eX_psi(:,:,k) / obj.eX_gamma(k)) / (obj.eX_gamma(k) - obj.dX - 1);
            end
            
            
            for n = 1: cols(Y)
                          
                for k = 1: obj.M
                    
                    %Auxiliary matrices
                    e_tt = latentStatistics.e_t(:,n,k)*latentStatistics.e_t(:,n,k)' + latentStatistics.cov_t * latentStatistics.cov_t_weights(n,k);
                    diagTraces = zeros(obj.dX,obj.dX);
                    
                    for j = 1 : obj.dY
                        diagTraces(j,j) = trace(e_tt * obj.covX_W(:,:,j,k));
                    end
                    
                    
                    %cluster specific contributions
                    eX_cov_pred(:, :, n) = eX_cov_pred(:, :, n) + latentStatistics.r(n,k) * (e_inv_u(n,k) * e_inv_psi(:,:,k) ...
                        + obj.eX_W(:,:,k) * latentStatistics.e_t(:,n,k) * obj.eX_mu(:,k)' ...
                        + obj.eX_mu(:,k) * latentStatistics.e_t(:,n,k)' * obj.eX_W(:,:,k)' ...
                        + obj.covX_mu(:,:,k) + obj.eX_mu(:,k)*obj.eX_mu(:,k)' ...
                        + obj.eX_W(:,:,k) *(e_tt) * obj.eX_W(:,:,k)' ...
                        + diagTraces ); 
                end
                %Covariances
                eX_cov_pred(:, :, n) = eX_cov_pred(:, :, n) - eX_mu_pred(:,n) * eX_mu_pred(:,n)';
            end
            
       end
        
        
        function [obj, energy] = updateParameterPosterior(obj,X,Y)
            %UPDATEPARAMETERPOSTERIOR - parameter posterior distribution
            %update routines, for internal use in learning
            
            % Auxiliary variables
            %------------------------------------------------------------
            
            %Scatter matrices
            tScatter = zeros(obj.D,obj.D,obj.M);
             
            
            %Kullback-Leibler divergences
            % Means
            DKL_muX = zeros(obj.M,1);
            DKL_muY = zeros(obj.M,1);
            
            %Projections
            DKL_WX = zeros(obj.M, obj.dX);
            DKL_WY = zeros(obj.M, obj.dY);
            
            %Precision
            DKL_psiX = zeros(obj.M,1);
            DKL_psiY = zeros(obj.M,1);
            
            % Means
            %------------------------------------------------------------
            for k = 1 : obj.M
                %Mean covariances
                obj.covX_mu(:,:,k) = inv( obj.betaX*eye(obj.dX) + ...
                    obj.eX_psi(:,:,k)* (obj.e_u(:,k)'* obj.r(:,k)));
                obj.covY_mu(:,:,k) = inv( obj.betaY*eye(obj.dY) + ...
                    obj.eY_psi(:,:,k) * (obj.e_u(:,k)' * obj.r(:,k)));
                
                %First moments
                obj.eX_mu(:,k) = obj.covX_mu(:,:,k)*(obj.betaX*obj.mX(:,k) + ...
                    obj.eX_psi(:,:,k) * row_sum(scale_cols((X - obj.eX_W(:,:,k)*obj.e_t(:,:,k)), ...
                    obj.e_u(:,k).*obj.r(:,k))));
                
                 obj.eY_mu(:,k) = obj.covY_mu(:,:,k)*(obj.betaY*obj.mY(:,k) + ...
                    obj.eY_psi(:,:,k) * row_sum(scale_cols((Y - obj.eY_W(:,:,k)*obj.e_t(:,:,k)), ...
                    obj.e_u(:,k).*obj.r(:,k))));
                
            end
            
            % Precision
            %------------------------------------------------------------
            
            % Degrees of freedom
            obj.eX_gamma = repmat(obj.gammaX,obj.M,1) + col_sum(obj.r)';
            obj.eY_gamma = repmat(obj.gammaY,obj.M,1) + col_sum(obj.r)';
            
            % The scale matrix
            for k = 1 : obj.M
                
                %Compute scatter matrices, also needed for the projection
                %matrix evaluations. 

                tScatter(:,:,k) = (repmat(sqrt(obj.r(:,k)'),obj.D,1) .* repmat(sqrt(obj.e_u(:,k)'),obj.D,1) .* ...
                    obj.e_t(:,:,k)) * (repmat(sqrt(obj.r(:,k)'),obj.D,1) .* repmat(sqrt(obj.e_u(:,k)'),obj.D,1) .* ...
                    obj.e_t(:,:,k))' + sum(obj.r(:,k))*obj.cov_t(:,:,k);
                          
                
                %Traces on the diagonals
                diagX = zeros(obj.dX,obj.dX);
                for j = 1 : obj.dX
                    diagX(j,j) = trace(tScatter(:,:,k) * obj.covX_W(:,:,j,k));
                end
                
                diagY = zeros(obj.dY,obj.dY);
                for j = 1 : obj.dY
                    diagY(j,j) = trace(tScatter(:,:,k) * obj.covY_W(:,:,j,k));
                end
                
                %Needed for other scatter matrices
                tmpX = repmat(sqrt(obj.r(:,k)'),obj.dX,1) .* repmat(sqrt(obj.e_u(:,k)'),obj.dX,1) .* ...
                    (X - obj.eX_W(:,:,k)*obj.e_t(:,:,k) - repmat(obj.eX_mu(:,k),1,obj.N));
            
                tmpY = repmat(sqrt(obj.r(:,k)'),obj.dY,1) .* repmat(sqrt(obj.e_u(:,k)'),obj.dY,1) .* ...
                    (Y - obj.eY_W(:,:,k)*obj.e_t(:,:,k) - repmat(obj.eY_mu(:,k),1,obj.N));
                
                
                %Moment evaluations ( inv -> /)
                obj.eX_psi(:,:,k) = obj.eX_gamma(k) * inv(obj.phiX + tmpX * tmpX' + obj.covX_mu(:,:,k)* ...
                    (obj.r(:,k)'*obj.e_u(:,k)) + sum(obj.r(:,k))* obj.eX_W(:,:,k) * obj.cov_t(:,:,k)* ...
                    obj.eX_W(:,:,k)' + diagX);
                
                obj.eY_psi(:,:,k) = obj.eY_gamma(k) * inv(obj.phiY + tmpY * tmpY' + obj.covY_mu(:,:,k)* ...
                    (obj.r(:,k)'*obj.e_u(:,k)) + sum(obj.r(:,k))* obj.eY_W(:,:,k) * obj.cov_t(:,:,k)* ...
                    obj.eY_W(:,:,k)' + diagY);
                    
            end
            
            
            % Projections
            %------------------------------------------------------------
            
            for k = 1 : obj.M
               
                for j = 1 : obj.dX
                    
                    %Covariances
                    obj.covX_W(:,:,j,k) = inv( diag(obj.e_alphaX(:,k)) * eye(obj.D) + obj.eX_psi(j,j,k)*tScatter(:,:,k));
                    
                    %Auxiliary matrix for the last part of the moment
                    tmpX = scale_rows(obj.eX_W(:,:,k),obj.eX_psi(j,:,k));
                    tmpX(j,:) = zeros(1,obj.D);
                    
                    %Means
                    obj.eX_W(j,:,k) = (obj.covX_W(:,:,j,k)*(diag(obj.e_alphaX(:,k)) * obj.VX(j,:,k)' + ...
                        row_sum(scale_cols(obj.e_t(:,:,k), obj.e_u(:,k).* obj.r(:,k) .* (obj.eX_psi(j,:,k)* ...
                        (X - repmat(obj.eX_mu(:,k),1,obj.N)))')) - tScatter(:,:,k) * row_sum(tmpX')))';
                    
                end
                
                for j = 1 : obj.dY
                    %Covariances
                    obj.covY_W(:,:,j,k) = inv( diag(obj.e_alphaY(:,k)) * eye(obj.D) + obj.eY_psi(j,j,k)*tScatter(:,:,k));
                    
                    %Auxiliary matrix for the last part of the moment
                    tmpY = scale_rows(obj.eY_W(:,:,k),obj.eY_psi(j,:,k));
                    tmpY(j,:) = zeros(1,obj.D);
                    
                    %Means
                    obj.eY_W(j,:,k) = (obj.covY_W(:,:,j,k)*(diag(obj.e_alphaY(:,k)) * obj.VY(j,:,k)' + ...
                        row_sum(scale_cols(obj.e_t(:,:,k), obj.e_u(:,k).* obj.r(:,k) .* (obj.eY_psi(j,:,k)* ...
                        (Y - repmat(obj.eY_mu(:,k),1,obj.N)))')) - tScatter(:,:,k) * row_sum(tmpY')))';
                    
                end

            end
            
            % ARD-priors
            %------------------------------------------------------------
                            
            obj.alphaX_a = repmat(obj.ard_a + obj.dX/2, obj.D, obj.M);
            obj.alphaY_a = repmat(obj.ard_a + obj.dY/2, obj.D, obj.M);

            for k = 1 : obj.M
                
                varWX = zeros(obj.dX, obj.D);
                varWY = zeros(obj.dY, obj.D);
                
                for j = 1 : obj.dX
                    varWX(j,:) = diag(obj.covX_W(:,:,j,k))';
                end
                
                for j = 1 : obj.dY
                    varWY(j,:) = diag(obj.covY_W(:,:,j,k))';
                end
                
                %TBD need to fix V variance if it is rv.
                
                obj.alphaX_b(:,k) = repmat(obj.ard_b,obj.D,1) + 0.5* ...
                    col_sum((obj.eX_W(:,:,k) - obj.VX(:,:,k)).*(obj.eX_W(:,:,k)-obj.VX(:,:,k)) + varWX)';
        
                obj.alphaY_b(:,k) = repmat(obj.ard_b,obj.D,1) + 0.5* ...
                    col_sum((obj.eY_W(:,:,k) - obj.VY(:,:,k)).*(obj.eY_W(:,:,k)-obj.VY(:,:,k)) + varWY)';
            end
            
            obj.e_alphaX = obj.alphaX_a ./ obj.alphaX_b;
            obj.e_alphaY = obj.alphaY_a ./ obj.alphaY_b;
            
            e_log_alphaX = digamma(obj.alphaX_a) - log(obj.alphaX_b);
            e_log_alphaY = digamma(obj.alphaX_a) - log(obj.alphaY_b);

            
            % Parameter KL-divergences
            %------------------------------------------------------------
            p_a = repmat(obj.ard_a, obj.D, obj.M);
            p_b = repmat(obj.ard_b, obj.D, obj.M);

            
            %ARD part KL-divergences
            DKL_alphaX = gammaln(p_a) - gammaln(obj.alphaX_a ) + obj.alphaX_a  .* log(obj.alphaX_b) - ...
                p_a .* log(p_b) + (obj.alphaX_a  - p_a).*(digamma(obj.alphaX_a ) - log(obj.alphaX_b)) + ...
                obj.alphaX_a  .* ((p_b - obj.alphaX_b) ./ obj.alphaX_b);
                    
            DKL_alphaY = gammaln(p_a) - gammaln(obj.alphaY_a ) + obj.alphaY_a  .* log(obj.alphaY_b) - ...
                p_a .* log(p_b) + (obj.alphaY_a  - p_a).*(digamma(obj.alphaY_a ) - log(obj.alphaY_b)) + ...
                obj.alphaY_a  .* ((p_b - obj.alphaY_b) ./ obj.alphaY_b);
            
            
            
            %Mean prior covariances
            covX_muprior = (1 / obj.betaX) * eye(obj.dX);
            covY_muprior = (1 / obj.betaY) * eye(obj.dY);
                        
            for k = 1 : obj.M
             
                % Means
                DKL_muX(k) = vbmcca.DKL_multinormal(obj.eX_mu(:,k),obj.covX_mu(:,:,k),obj.mX(:,k),covX_muprior);
                DKL_muY(k) = vbmcca.DKL_multinormal(obj.eY_mu(:,k),obj.covY_mu(:,:,k),obj.mY(:,k),covY_muprior);
                
                %Projections (additional terms are due to gamma logarithm expectation)
                for j = 1 : obj.dX
                    DKL_WX(k,j) = vbmcca.DKL_multinormal(obj.eX_W(j,:,k)', obj.covX_W(:,:,j,k), ...
                        obj.VX(j,:,k)', inv(diag(obj.e_alphaX(:,k)))) ...
                        - 0.5*sum(e_log_alphaX(:,k)) + 0.5*sum(log(obj.e_alphaX(:,k)));
                end
                
                for j = 1 : obj.dY
                    DKL_WY(k,j) = vbmcca.DKL_multinormal(obj.eY_W(j,:,k)', obj.covY_W(:,:,j,k), ...
                        obj.VY(j,:,k)', inv(diag(obj.e_alphaY(:,k)))) ...
                        - 0.5*sum(e_log_alphaY(:,k)) + 0.5*sum(log(obj.e_alphaY(:,k)));
                end
                    
                %Precision
                DKL_psiX(k) = vbmcca.DKL_wishart(obj.eX_psi(:,:,k), obj.eX_gamma(k), obj.gammaX  * inv(obj.phiX), ...
                    obj.gammaX,obj.dX);
                
                DKL_psiY(k) = vbmcca.DKL_wishart(obj.eY_psi(:,:,k), obj.eY_gamma(k), obj.gammaY  * inv(obj.phiY), ...
                    obj.gammaY,obj.dY);
            
            end
            energy = (-1) * (sum(DKL_muX + DKL_muY + DKL_psiX + DKL_psiY + ...
                row_sum(DKL_WX) + row_sum(DKL_WY)) + sum(row_sum(DKL_alphaY) + row_sum(DKL_alphaX)));
        end
        
            
        
        function [model, lowerbound] = learn(obj, X, Y, parameters)
        %LEARN - VBEM iteration; Initialise learning data structures first
        %with for isntance kmeans initialisation
        %
        % learn(obj, X, Y, parameters) - parameters:
        % X - data matrix X, dim(X) = dX X N
        % Y - data matrix Y, dim(Y) = dY X N
        % parameters - learningParameters struct 
        %
        % output:
        % model - the trained model
        % lowerbound - a table of the variational lowerbounds
        
            penergies = zeros(parameters.maxIter, 1);
            lenergies = zeros(parameters.maxIter, 1);
            
            %Sweep latents once after k-means init to get all expectations
            %evaluated, 0 to nu-updates ; 1 to regularised cluster
            %determination
            
            [model,foo] = obj.updateLatentPosterior(X,Y,0,1);
            
            %Iterate until maxIter or tolerance boundary is reached
            wnew = [];
            for i = 1:parameters.maxIter
              
                if i>1
                    wold = wnew;
                end

                [model,penergies(i)] = model.updateParameterPosterior(X, Y);
                [model,lenergies(i)] = model.updateLatentPosterior(X, Y, parameters.nuUpdates, i < 10);
              
                % Edit by Will Penny, to print out progress only every 16
                % iterations
                if mod(i,16)==0
                    disp(strcat('Iteration ', num2str(i),', Lower bound: ', num2str(penergies(i)+lenergies(i))));
                end
                
                %Convergence test
                % if( i >= 2 && abs( (penergies(i)+lenergies(i)) - (penergies(i-1)+lenergies(i-1)) ) ...
                %         / abs(penergies(i)+lenergies(i)) < parameters.relativeTol)
                % 
                %     disp(strcat('Relative change smaller than the tolerance after ',num2str(i),' iterations.'));
                %     break;
                % end

                wnew = [model.eX_W(:);model.eY_W(:)];
                if i > 8 
                    changeW=norm(wnew-wold)/length(wold);
                    if changeW < parameters.relativeTol
                        disp(strcat('Relative change smaller than the tolerance after ',num2str(i),' iterations.'));
                        break;
                    end
                end
            end
            lowerbound = penergies(1:i) + lenergies(1:i);
        end
        
        function [Ux,Uy,cor] = findRotation(obj,k)
        %FINDROTATION - Returns the nonprobabilistic CCA-solution by
        %solving the rotatiol ambiguity (see "Robust probabilistic projections" paper for details)
        %
        % [Ux,Uy,cor] = findRotation(obj,k)
        %
        % parameters:
        % k - the cluster index from which the parameters are calculated
        %
        % output:
        % Ux - a Matrix where the columns are CCA-projections for dataset X
        % Uy - a Matrix where the columns are CCA-projections for dataset Y
        % cor - diagonal matrix of the estimated canonical correlations
        
            Bx = obj.eX_W(:,:,k)' * obj.eX_psi(:,:,k) * obj.eX_W(:,:,k) + eye(obj.D);
            By = obj.eY_W(:,:,k)' * obj.eY_psi(:,:,k) * obj.eY_W(:,:,k) + eye(obj.D);

            Rt1 = sqrtm(eye(obj.D) - inv(Bx)) * (eye(obj.D)-inv(By)) * sqrtm(eye(obj.D)-inv(Bx));
            Rt2 = sqrtm(eye(obj.D) - inv(By)) * (eye(obj.D)-inv(Bx)) * sqrtm(eye(obj.D)-inv(By));

            [R1,D1] = eig(Rt1);
            [R2,D2] = eig(Rt2);
  
            covX = obj.eX_W(:,:,k)*obj.eX_W(:,:,k)' + inv(obj.eX_psi(:,:,k));
            covY = obj.eY_W(:,:,k)*obj.eY_W(:,:,k)' + inv(obj.eY_psi(:,:,k));

            isqX = inv(sqrtm(eye(obj.D) - inv(Bx)));
            isqY = inv(sqrtm(eye(obj.D) - inv(By)));

            Ux = covX \ obj.eX_W(:,:,k) * isqX * R1;
            Uy = covY \ obj.eY_W(:,:,k) * isqY * R2;
            
            cor = sqrtm(D1);
            
            %Sort wrt correlations
            [cor, perm] = sort(diag(cor),'descend');
            Ux = Ux(:,perm);
            Uy = Uy(:,perm);
        end 
        
        function statistics = defaultComparisonStatistics(obj)
        % DEFAULTCOMPARISONSTATISTICS - Returns the learned parameters
        % as the comparison statistics
        %
        % statistics = defaultComparisonStatistics(obj)
        % 
        % output:
        % statistics - comparisonStatistics structure
            statistics = obj.initLatentStatistics(obj.D, obj.N, obj.M);
            statistics.e_t = obj.e_t;
            statistics.cov_t = obj.cov_t;
            statistics.cov_t_weights = 1 ./ obj.e_u;
            statistics.r = obj.r;
        end
    end
    
    % Auxiliary static methods
    %
    % Divergence methods use some shortcuts based on notations/results in
    % derivations so they might not be usable for general stuff without
    % minor modifications
    
    methods(Static)
        
        function parameters = initLearningParameters(maxIter, tol)
        %INITLEARNINGPARAMETERS - Creates a structure of learning parameters for
        %variational learning method
        %
        % parameters = initLearningParameters(maxIter, tol)
        %
        % parameters:
        % maxIter - Maximum number of the VBEM iterations
        % tol - Iteration is terminated if the relative change in the
        % log-likelihoods is smaller than tol.
        % nuUpdates - Whether hyperparameter nu is updated (1) or not (0).
        %
        % output:
        % parameters - a structure that encapsulates the learning settings
        
            parameters.maxIter = maxIter;
            parameters.relativeTol = tol;
            parameters.nuUpdates = 1;
            %TBD additional knobs
        end
        
        function statistics = initLatentStatistics(D, N, M)
        % INITLATENSTATISTICS - A method for setting up a structure
        % which consists of sufficient latent statistics. These statistics can be
        % calculated by, for instance, methods inferLatents(X/Y). The
        % end used should not have any need to initialise the struct
        % manually.
        
            statistics.N = N;
            statistics.M = M;
            statistics.D = D;
            statistics.e_t = zeros(D,N,M);
            statistics.u_alpha = zeros(N, M);
            statistics.u_beta =  zeros(N, M);
            statistics.cov_t = zeros(D,D,M);
            statistics.cov_t_weights = ones(N,M);
            statistics.r = zeros(N,M);
        end
        
        function div = DKL_multinormal(mu0,cov0,mu1,cov1)
        %DKL_MULTINORMAL - Kullback-Leibler divergence for multivariate
        %normal distribution DKL(\mathcal{N0},\mathcal{N1})
        
            div =  0.5 * (logdet( cov1 ) - logdet( cov0 ) + ...
                    trace(inv(cov1)*cov0) + ...
                    ( mu1 - mu0 )'*inv(cov1)*( mu1 - mu0) - length(mu1)); 
               
        end
        
        function val = eLogDetPsi(ePsi,p,D)
        %ELOGDETPSI - Expectation for logarithm of determinant of wishart
        %distributed random matrix, ePsi is expectation of corresponding
        %wishart rv.
        
        val = sum(digamma((p + 1 -(1:D))/2)) + ...
                    D*log(2) + logdet(ePsi * 1/p);
        end
       
        function val = logGenGamma(p,D)
        %LOGGENGAMMA - logarithm of the generalised gamma function,
        %WITHOUT(!) normalisation constant
        
            val = sum(gammaln((p+1-(1:D))/2)); 
        
        end
        
        function div = DKL_wishart(ePsi0,gamma0,ePsi1,gamma1,D)
        %DKL_WISHART - Kullback-leibler divergence for Wishart distribution
        
            div = (gamma0-1 -D)/2 * vbmcca.eLogDetPsi(ePsi0,gamma0,D) - ...
                  (gamma1-1 -D)/2 * vbmcca.eLogDetPsi(ePsi1,gamma1,D)- gamma0*D/2 + ...
                  gamma0/2 * trace(inv(1/gamma1 * ePsi1) * ePsi0 * 1/gamma0) + ...
                  (gamma1-gamma0)*D/2*log(2) + vbmcca.logGenGamma(gamma1,D) - vbmcca.logGenGamma(gamma0,D) + ...
                  gamma1/2*logdet(ePsi1 * 1/gamma1) - gamma0/2*logdet(ePsi0 * 1/gamma0);
                  
        end     
    end
end

