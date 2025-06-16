function [] = blr_plot_boundary (bsr,x)
% Plot decision boundaries for 2D data
% FORMAT [] = blr_plot_boundary (blr,x)
%
% M         BLR data structure
% x         [d x N] input data with N exemplars

[N,d]=size(x);
if d > N, x=x'; end
xmin=min(x(:,1));
xmax=max(x(:,1));
dx=(xmax-xmin)/100;
u=[xmin:dx:xmax];

w = bsr.m;
m = -w(1)/w(2);
c = (0.5-w(3))/w(2);
y = m*u+c;
hold on
plot(u,y,'k-');
