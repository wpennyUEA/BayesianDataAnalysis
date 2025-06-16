function [] = mvn_plot2D (m,C,R)
% Plot 2D multivariate normal density
% FORMAT [] = mvn_plot2D (m,C,R)
%
% m,C           Parameters of multivariate normal
% R             Definition of plotting region

sx=diag(sqrt(C));
nsd=4;

try N1=R.N1; catch N1=100; end
try N2=R.N2; catch N2=N1; end
try x1_min=R.x1_min; catch x1_min=m(1)-nsd*sx(1); end
try x1_max=R.x1_max; catch x1_max=m(1)+nsd*sx(1); end
try x2_min=R.x2_min; catch x2_min=m(2)-nsd*sx(2); end
try x2_max=R.x2_max; catch x2_max=m(2)+nsd*sx(2); end

x1=linspace(x1_min,x1_max,N1);
x2=linspace(x2_min,x2_max,N2);
x=[];
for j=1:N2,
    for i=1:N1,
        xt=[x1(i),x2(j)]';
        x=[x,xt];
    end
end
p = spm_mvNpdf (x,m,C);
p = reshape(p,N1,N2);
imagesc(x1,x2,log(p'),[-20 0]);
axis xy
set(gca,'FontSize',16);
xlabel('x_1');
ylabel('x_2');