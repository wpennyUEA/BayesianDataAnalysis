function [h,p] = vbcca_query (cca,org,scaling,con)
% Compute surprise of query data under CCA normative model
% FORMAT [h,p] = vbcca_query (cca,org,scaling,con)
%
% cca           normative model
% org           .X1 and .X2 data in original space 
% scaling       .msd1 and msd2 means and SDs of original data
% con           contrast matrices (optional)
%               .Gamma1 and .Gamma2 (see vbcca_cond_subspace)
% 
% h             surprise
% p             probability

X1 = norm_space (org.X1,scaling.msd1);
X2 = norm_space (org.X2,scaling.msd2);

if nargin>3
    subspace = 1;
    c1 = con.Gamma1*X1;
    c2 = con.Gamma2*X2;
else
    subspace = 0;
end

Ntest = size(X1,2);
for n=1:Ntest,
    if subspace        
        [tmp1,tmp2,p(n)] = vbcca_cond_subspace(cca,c2(:,n),con,c1(:,n));
    else
        p(n) = vbcca_conditional (cca,X1(:,n),X2(:,n));
    end
end
h = -log(p);
