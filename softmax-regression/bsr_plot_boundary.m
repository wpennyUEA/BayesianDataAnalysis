function [] = bsr_plot_boundary (bsr,x)
% Plot decision boundaries for 2D data
% FORMAT [] = bsr_plot_boundary (bsr,x)
%
% bsr       bsr data structure
% x         [d x N] input data with N exemplars

[N,d]=size(x);
if d > N, x=x'; end
xmin=min(x(:,1));
xmax=max(x(:,1));
dx=(xmax-xmin)/100;
u=[xmin:dx:xmax];

K=length(bsr.class);
for k1=1:K,
    w1=bsr.class(k1).w;
    for k2=k1:K,
        % plot decision boundaries between categories k1 and k2
        w2=bsr.class(k2).w;
        denom = w1(2)-w2(2);
        m = (-w1(1)+w2(1))/denom;
        c = (-w1(3)+w2(3))/denom;
        y = m*u+c;
        hold on
        plot(u,y,'r-');
    end
end
