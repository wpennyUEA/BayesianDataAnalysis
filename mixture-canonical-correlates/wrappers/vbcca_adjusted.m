function [X1adj,X2adj,c1hat] = vbcca_adjusted (cca,con,X1,X2)
% Adjust data for contrasts
% FORMAT [X1adj,X2adj,c1hat] = vbcca_adjusted (cca,con,X1,X2)
%
% cca       fitted model
% con       .Gamma1 and .Gamma2 are contrast matrices
% X1        original X1 data
% X2        original X2 data
%
% X1adj     adjusted X1 data
% X2adj     adjusted X2 data
% c1hat     the adjustment for X1

X2adj = con.Gamma2*X2;

d2 = length(con.Gamma2);
loc1 = find(con.Gamma2==1);
ind = [1:d2];
ind(loc1)=[];
if length(loc1)>1 | length(loc1)==0 | any(~(con.Gamma2(ind)==0))
    disp('Error in vbcca_adjusted.m: con.Gamma2 must be a zero vector with a single 1 entry');
    disp('con.Gamma2 is:');
    disp(con.Gamma2);
    return
end

% Compute R2 - the complement of Gamma2
R2 = eye(d2);
R2(loc1,:)=[];

N = size(X2,2);
for n=1:N,
    c2 = R2*X2(:,n);
    c1hat(n) = vbcca_cond_subspace (cca,c2,con);
end

% Subtract from X1 data that which can be predicted from other X2 variables
X1adj = con.Gamma1*X1 - c1hat;
%X1adj = con.Gamma1*X1;


