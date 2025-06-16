function [] = mvt_plot2D (mu,Lambda,v,R)
% Plot 2D multivariate T density
% FORMAT [] = mvt_plot2D (mu,Lambda,v,R)
%
% mu,Lambda,v   Parameters of multivariate T
% R             Definition of plotting region
%__________________________________________________________________________
% Copyright (C) 2014 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id$

Cx=(v/(v-2))*inv(Lambda);
sx=diag(sqrt(Cx));
nsd=4;

% if nargin < 4 | isempty(R)
%     R=[];
% end

try N1=R.N1; catch N1=100; end
try N2=R.N2; catch N2=N1; end
try x1_min=R.x1_min; catch x1_min=mu(1)-nsd*sx(1); end
try x1_max=R.x1_max; catch x1_max=mu(1)+nsd*sx(1); end
try x2_min=R.x2_min; catch x2_min=mu(2)-nsd*sx(2); end
try x2_max=R.x2_max; catch x2_max=mu(2)+nsd*sx(2); end

x1=linspace(x1_min,x1_max,N1);
x2=linspace(x2_min,x2_max,N2);
x=[];
for j=1:N2,
    for i=1:N1,
        xt=[x1(i),x2(j)]';
        x=[x,xt];
    end
end
p = spm_mvtpdf (x,mu,Lambda,v);
p = reshape(p,N1,N2);
imagesc(x1,x2,log(p'),[-20 0]);
axis xy
set(gca,'FontSize',16);
xlabel('x_1');
ylabel('x_2');
