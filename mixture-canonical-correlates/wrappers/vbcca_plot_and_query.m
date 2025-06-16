function [R] = vbcca_plot_and_query (cca,con,Q,P)
% Plot normative density and query data and compute surprise thereof
% FORMAT [R] = vbcca_plot_and_query (cca,con,Q,P)
%
% ------------------------------------------------------------------------
% INPUTS:
% cca           normative model
%
% con           contrast matrices 
%                  .Gamma1 and .Gamma2 (see vbcca_cond_subspace)
%
% Q             Query Info
%                   .X1 and .X2 data    in original space 
%                   .msd1 and msd2      means and SDs of original reference data
%                   .percentile         threshold for anomalies
%                   .c2_delta           c2 variable must be within this tolerance
%                                       to be plotted
%
% P             Plotting Info
%                   .c1_range           range(s) of contrasted state variable(s)
%                   .c2                 contrasted context variable
%                   .c1_name            name(s) of contrasted state variable(s)
%                   .c2_name            name of context variable
%                   .c1_offset          offset state values in plots
%                   .plot_mode = 1 as default
%                   .plot_mean = 1 as default
%
% ------------------------------------------------------------------------
% OUTPUTS:
% R             Report Info
%                   .h             surprise
%                   .p             probability
%                   .qind          indices of anomalous data points 
%                                  (with h in a percentile above
%                                  Q.percentile)
%                   .outscore      univariate outlier scores
%
% If Q supplied, compute surprise of query data
% If P supplied, plot normative density
% If Q and P supplied, add anomalous query data points to plot

try plot_mode = P.plot_mode; catch plot_mode=1; end
try plot_mean = P.plot_mean; catch plot_mean=1; end

%----------------------------------------------------------------------
% Compute surprise of query data
X1 = norm_space (Q.X1,Q.msd1);
X2 = norm_space (Q.X2,Q.msd2);

subspace = 1;
c1 = con.Gamma1*X1;
c2 = con.Gamma2*X2;
n1 = size(con.Gamma1);
query_c1_org = con.Gamma1*Q.X1;
query_c2_org = con.Gamma2*Q.X2;

Ntest = size(X1,2);
for n=1:Ntest,
    if subspace
        [x1_hat,gamma,p(n),T,std_dev] = vbcca_cond_subspace(cca,c2(:,n),con,c1(:,n));
        if length(cca.pi)==1
            % single cluster model
            outscore(:,n) = (c1(:,n)-x1_hat)./std_dev;
        else
            K = length(cca.pi);
            for k=1:K,
                outsc(:,k) = (c1(:,n)-x1_hat)./std_dev(:,k);
            end
            outscore(:,n) = outsc*gamma;
        end
    else
        p(n) = vbcca_conditional (cca,X1(:,n),X2(:,n));
    end
end
R.p = p;
R.h = -log(p);
R.qind = [];
R.outscore = outscore;

if nargin < 4
    return
end

[tmp,ind] = sort(R.h);
k = ceil(Ntest*(Q.percentile/100));
k = max([1,k]);
qind = ind(k:Ntest);
R.qind = qind;

%----------------------------------------------------------------------
% Plot Normative Density

msd1 = con.Gamma1*Q.msd1;
msd2 = con.Gamma2*Q.msd2;
n1 = size(con.Gamma1,1);

switch n1
    case 1,
        % 1D dependent variable
        bins = 100;
        norm_range = norm_space(P.c1_range,msd1);
        S.xmin = norm_range(1,1);
        S.xmax = norm_range(1,2);
        S.dx = (S.xmax-S.xmin)/bins;
        c1 = [S.xmin:S.dx:S.xmax]';
        c1_org = org_space(c1(:)',msd1);
        
        norm_c2 = norm_space (P.c2,msd2);
        N = length(norm_c2);
        if size(con.Gamma1,1) == 1
            for n=1:N,
                [c1_pred(n),tmp2,p1(n,:)] = vbcca_cond_subspace (cca,norm_c2(n),con,c1);
            end
        else
            disp('Error in vbcca_plot_and_query.m: x1 subspace expected to be 1D');
            return
        end
        
        figure
        imagesc(P.c2,c1_org+P.c1_offset,p1');
        axis xy
        colormap gray
        xlabel(P.c2_name);
        ylabel(P.c1_name);
        colormap gray
        hold on
        
        if plot_mode
            for n=1:N,
                [tmp,ind]=max(p1(n,:));
                c1_mode(n) = c1_org(ind);
            end
            plot (P.c2,c1_mode+P.c1_offset,'b');
        end
        
        if plot_mean
            c1_org = org_space(c1_pred,msd1);
            plot(P.c2,c1_org+P.c1_offset,'r');
        end
        
        
        %----------------------------------------------------------------------
        % Also plot query data points above percentile threshold
        hold on
        plot(query_c2_org(qind),query_c1_org(qind)+P.c1_offset,'rx');
        
        
    case 2,
        % 2D dependent variable
        bins = 30;
        norm_range = norm_space (P.c1_range,msd1);
        S.xmin = norm_range(2,1);
        S.xmax = norm_range(2,2);
        S.ymin = norm_range(1,1);
        S.ymax = norm_range(1,2);
        S.dx = (S.xmax-S.xmin)/bins;
        S.dy = (S.ymax-S.ymin)/bins;
        
        x = [S.xmin:S.dx:S.xmax];
        y = [S.ymin:S.dy:S.ymax];
        
        org = org_space([x;y],msd1);
        
        Nx = length(x);
        Ny = length(y);
        
        k = 1;
        for i=1:Nx,
            for j=1:Ny,
                X1p(:,k)= [x(i) y(j)]';
                k = k+1;
            end
        end
        
        norm_c2 = norm_space (P.c2,msd2);
       
        Nc2 = length(norm_c2);
        figure
        rNc2 = ceil(sqrt(Nc2));
        for m = 1:Nc2,
            [tmp1,tmp2,p1cond] = vbcca_cond_subspace (cca,norm_c2(m),con,X1p);
            
            k = 1;
            for i=1:Nx,
                for j=1:Ny,
                    D1(j,i)=p1cond(k);
                    k = k+1;
                end
            end
            
            subplot(rNc2,rNc2,m);
            imagesc(org(1,:)+P.c1_offset(1),org(2,:)+P.c1_offset(2),D1);
            xlabel(P.c1_name{1});
            ylabel(P.c1_name{2});
            colormap gray
            axis xy
            hold on
            title(sprintf('%s = %1.2f',P.c2_name,P.c2(m)));
            
            %----------------------------------------------------------------------
            % Plot anomalous query data points within Q.c2_delta of P.c2(1)          
            a = query_c1_org(1,qind);
            b = query_c1_org(2,qind);
            ind = find(abs(query_c2_org(qind)-P.c2(m)) < Q.c2_delta);
            plot(a(ind)+P.c1_offset(1),b(ind)+P.c1_offset(2),'rx');
        end
        
end