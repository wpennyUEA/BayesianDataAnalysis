
close all
clear all

% Load Image
x=imread('Alan','jpg');

% Crop Alan
startx=100;
starty=300;

N=256;
xg=double(x(startx:startx+N-1,starty:starty+N-1,2));
xg=xg-mean(mean(xg));

% Show Alan
figure
imagesc(xg);
colormap gray
axis image
title('I''m Alan Partridge');

disp('Doing SVD .... !');
[u,s,v]=svd(xg,0);
eval=diag(s);
var=eval.^2;
var_exp=cumsum(var)/sum(var);
figure
plot(var_exp);
title('Variance Explained');
grid on
xlabel('Number of Components');

% Reconstruct 
comps=20;

figure
cyber_alan=zeros(N,N);
for i=1:comps,
    cyber_alan=cyber_alan+eval(i)*u(:,i)*v(:,i)';
    imagesc(cyber_alan);
    colormap gray
    axis image
    title(sprintf('With %d components',i));
    disp('Press a key ...');
    pause;
end

title(sprintf('Aha ! I''m Cyber-Alan ! Created with %d singular components',comps));
