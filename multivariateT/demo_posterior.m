
clear all
close all

% Normal-Wishart model parameters
P=2;
a0=P/2;
B0=diag([1 1]);
beta0=0.01;
m0=zeros(P,1);

% Plot prior predictive density
mu_w0=m0;
w_s=(beta0/(beta0+1))*(a0-0.5*(P-1));
Lambda_w0=w_s*inv(B0);
v_w0=2*a0-P+1;
figure
subplot(1,2,1);
mvt_plot2D (mu_w0,Lambda_w0,v_w0);
axis square
title('n=0');
pause(0.1);

M.prior.P=P;
M.prior.a=a0;
M.prior.B=B0;
M.prior.beta=beta0;
M.prior.m=m0;

% Generate data with known mu, C
new_data=0;
if new_data
    mu=[10,7]';
    s1=2;
    s2=0.5;
    r=-0.7;
    %r=0;
    c12=r*s1*s2;
    C=[s1^2 c12; c12 s2^2];
    Lambda=inv(C);
    
    N=32;
    x = spm_normrnd(mu, C, N);
    save xdata x N s1 s2 r mu C
else
    load xdata
end

R.x1_min=-5;
R.x1_max=20;
R.x2_min=3;
R.x2_max=10;

% Plot evolving posterior predictive density
for n=1:N,
    M = spm_nwpost (M,x(:,1:n));
    clf;
    
    subplot(1,2,1);
    mvt_plot2D (M.post.mu_w,M.post.Lambda_w,M.post.v_w,R);
    axis square
    hold on
    for j=1:n,
        plot(x(1,j),x(2,j),'kx','MarkerSize',10);
    end
    title(sprintf('MV-T, n=%d',n));
    
    if n>2
        sx=x(:,1:n);
        mw=mean(sx,2);
        Sw=cov(sx',1);
        subplot(1,2,2);
        mvn_plot2D (mw,Sw,R);
        axis square
        hold on
        for j=1:n,
            plot(x(1,j),x(2,j),'kx','MarkerSize',10);
        end
        title(sprintf('MV-Gauss, n=%d',n));
    end
    
%     disp('True Cov:');
%     C
%     if n>2
%     disp('Sample Cov:');
%     Sw
%     end
    disp('Posterior Mean NW-Cov:');
    M.post.B/M.post.a
    
    drawnow
    %pause(0.1);
end

% Maximum-Likelihood Estimates
sml_1=sqrt(Sw(1,1));
sml_2=sqrt(Sw(2,2));
rml=Sw(1,2)/(sml_1*sml_2);
    
% Sample precisions from Wishart prior
Ns=1000;
L=spm_wishrnd(M.prior.B,M.prior.a,Ns);
for s=1:Ns,
    C=inv(squeeze(L(:,:,s)));
    sig_1(s)=sqrt(C(1,1));
    sig_2(s)=sqrt(C(2,2));
    rw(s)=C(1,2)/(sig_1(s)*sig_2(s));
end

figure
plot(sig_1,sig_2,'.');
hold on
set(gca,'FontSize',18);
xlabel('\sigma_1');
ylabel('\sigma_2');
title('Prior');


% Sample precisions from Wishart posterior
Ns=1000;
L=spm_wishrnd(M.post.B,M.post.a,Ns);
for s=1:Ns,
    C=inv(squeeze(L(:,:,s)));
    sig_1(s)=sqrt(C(1,1));
    sig_2(s)=sqrt(C(2,2));
    rw(s)=C(1,2)/(sig_1(s)*sig_2(s));
end

figure
plot(sig_1,sig_2,'.');
hold on
grid on
plot(s1,s2,'rx','MarkerSize',20,'LineWidth',2);
plot(sml_1,sml_2,'gx','MarkerSize',20,'LineWidth',2);
set(gca,'FontSize',18);
xlabel('\sigma_1');
ylabel('\sigma_2');
title('Posterior');

figure
[n,c]=hist(rw,20);
n=n/sum(n);
bar(c,n);
mn=max(n);
set(gca,'FontSize',18);
xlabel('r');
hold on
plot([rml rml],[0 mn],'g','LineWidth',4);
plot([r r],[0 mn],'r','LineWidth',4);
grid on
ylabel('p(r|w)');