
close all
clear all

% This script uses Singular Value Decomposition (SVD) for image compression
% and dimensionality reduction

% Load Image
x=imread('Alan','jpg');

% Crop
startx=100;
starty=300;
N=256;  % Size of cropped square image (pixels)

% Extract green channel and convert to double precision
xg=double(x(startx:startx+N-1,starty:starty+N-1,2));
% Centre the image data
xg=xg-mean(mean(xg));

% Show cropped, mean-centred image
figure
imagesc(xg);
colormap gray
axis image
title('I''m Alan Partridge');
disp('Doing SVD .... !');

% Perform SVD on the image matrix
[u,s,v]=svd(xg,0);  % u = left singular vectors, s = diagonal matrix with singular vectors, v = right singular vectors
% Extract singular values
eval=diag(s);
% Calculate variance explained by each singular value
var=eval.^2;
% Calculate cumulative variance explained
var_exp=cumsum(var)/sum(var);

% Plot cumulative variance explained vs number of components
figure
plot(var_exp);
title('Variance Explained');
grid on
xlabel('Number of Components');

% Number of components in reconstruction
comps=20; 

% Reconstruct image incrementally
figure
cyber_alan=zeros(N,N);
for i=1:comps
    cyber_alan=cyber_alan+eval(i)*u(:,i)*v(:,i)';
    % Display current reconstruction
    imagesc(cyber_alan);
    colormap gray
    axis image
    title(sprintf('With %d components',i));
    disp('Press a key ...');
    pause;
end

% Final title
title(sprintf('Aha ! I''m Cyber-Alan ! Created with %d singular components',comps));
