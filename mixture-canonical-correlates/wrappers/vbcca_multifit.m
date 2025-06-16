function [M,DZ] = vbcca_multifit (Kmax,Mmax,D,verbose)
% Fit multiple VBCCA models to data
% FORMAT [M,DZ] = vbcca_multifit (Kmax,Mmax,D,verbose)
%
% Kmax      maximum number of latent vars
% Mmax      maximum number of clusters
% D         data to be fitted 
%           .X1, .X2, .X1names, .X2names
% verbose   (1) to plot results
%
% This function searches over 1..Kmax, 1..Mmax for model
% with the highest evidence
%
% M         Output data structure
%           .cca            best model
%           .res_cca(m,k)   model with m clusters and k latent vars
%           .F(m,k)         model evidences
%           .max_m          number of clusters in best model
%           .max_k          number of latents in best model
%           .logBF          log Bayes Factor of best model versus 
%                           best single cluster model
%
% DZ        normalised data 
%           .X1, .X2 
%           .m1,.m2,.s1,.s1 means and SDs of original data

try verbose = verbose; catch verbose=0; end

% Normalise all variables to zero mean, unit variance
[X1,DZ.m1,DZ.s1] = zmuv(D.X1);
[X2,DZ.m2,DZ.s2] = zmuv(D.X2);
DZ.X1 = X1;
DZ.X2 = X2;
DZ.X1names = D.X1names;
DZ.X2names = D.X2names;

Fmax = -Inf;
for m = 1:Mmax,
    for k=1:Kmax,
        cca = vbcca (X1,X2,k,m);
        res_cca(m,k) = cca;
        F(m,k) = cca.F;
        if cca.F > Fmax
            Fmax = cca.F;
            max_m = m;
            max_k = k;
        end
    end
end

% The best model
cca = res_cca(max_m,max_k);

M.res_cca = res_cca;
M.cca = cca;
M.F = F;
M.Fmax = Fmax;
M.max_m = max_m;
M.max_k = max_k;
M.logBF = Fmax - max(F(1,:)); % best overall vs. the best single cluster model

if verbose
    p={'One-Cluster','Two-Clusters','Three Clusters','Four Clusters'};
    col={'b-','k-','r-','g-'};
    for m=1:min([Mmax,4]),
        plot(F(m,:),col{m});
        hold on
    end
    xlabel('Latent Dimension');
    legend(p{1:Mmax});
end
