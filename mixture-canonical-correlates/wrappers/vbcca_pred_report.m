function [] = vbcca_pred_report (cca,DZ)
% Report on predictability of individual variables
% FORMAT [] = vbcca_pred_report (cca,DZ)
%
% cca       returned by vbcca.m
% DZ        normalised data (output of vbcca_multifit.m)

X1 = DZ.X1;
X2 = DZ.X2;
d1 = size(X1,1);
d2 = size(X2,1);

con.Gamma2 = eye(d2);
for d = 1:d1,
    con.Gamma1=zeros(1,d1);
    con.Gamma1(d)=1;
    
    c2 = con.Gamma2*X2;
    
    empirical = con.Gamma1*X1;
    [fitted,tmp1,tmp2,T] = vbcca_cond_subspace (cca,c2,con);
    dep(d).T = T;
    [R,P] = corrcoef(empirical(:),fitted(:));
    r(d) = R(1,2);
    pval(d) = P(1,2);
end

for d=1:d1,
    disp(sprintf('%s r = %1.3f, p-val = %1.3f',DZ.X1names{d},r(d),pval(d)));
end

pthresh = 0.1;
ind = find(pval<pthresh);
if length(ind) == 0
    disp(sprintf('No X1 variables predicted by X2 variables at p=%1.2f significance level',pthresh));    
    return
end


for i = 1:length(ind),
    d = ind(i);
    con.Gamma1=zeros(1,d1);
    con.Gamma1(d)=1;
    empirical = con.Gamma1*X1;
    
    % Plot univariate mappings
    h = figure;
    set(h,'Name',DZ.X1names{d});
    for j = 1:d2,
        con.Gamma2=zeros(1,d2);
        con.Gamma2(j)=1;
        c2 = con.Gamma2*X2;
        [c2s,sind] = sort(c2);
        fitted = vbcca_cond_subspace (cca,c2s,con);
        subplot(d2,1,j);
        plot(c2s,fitted,'k-');
        hold on
        plot(c2s,empirical(sind),'x');
        xlabel(DZ.X2names{j});
        ylabel(DZ.X1names{d});
        grid on
    end
    
    % Plot transformation matrices
    disp(' ');
    disp('Cluster Frequencies and Transformation matrices');
    for m=1:cca.M
        disp(sprintf('Cluster %d, freq = %1.2f, T:',m,cca.pi(m)));
        T = dep(d).T{m};
        disp(T);
    end
end