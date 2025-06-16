function [h2] = vbcca_marginal_2D (cca,S,X1,X2)
% Plot 2D marginal density 
% FORMAT [h2] = vbcca_marginal_2D (cca,S,X1,X2)
%
% cca       cca data structure from vbcca.m
% S         space to plot over
%           .xmin, .xmax, .dx
%           .ymin, .ymax, .dy
% X1        optionally overlay X1 data points 
% X2        optionally overlay X2 data points 
%           (if these fields included)
% 
% h2        figure pointer to plot of p(x2)


x = [S.xmin:S.dx:S.xmax];
y = [S.ymin:S.dy:S.ymax];
Nx = length(x);
Ny = length(y);

k = 1;
for i=1:Nx,
    for j=1:Ny,
        Ypred(:,k)= [x(i) y(j)]';
        k = k+1;
    end
end

X1p = Ypred;
X2p = Ypred;

[p1,p2] = vbcca_marginal (cca,X1p,X2p);

k = 1;
for i=1:Nx,
    for j=1:Ny,
        D1(j,i)=p1(k);
        D2(j,i)=p2(k);
        k = k+1;
    end
end

figure
imagesc(x,y,D1);
axis xy
if nargin>2
    hold on
    plot(X1(1,:),X1(2,:),'wo');
end
title('p(x1)');
xlabel('x1[1]');
ylabel('x1[2]');

figure;
imagesc(x,y,D2);
axis xy
if nargin>2
    hold on
    plot(X2(1,:),X2(2,:),'wo');
end
title('p(x2)');
xlabel('x2[1]');
ylabel('x2[2]');
h2 = gca;

 