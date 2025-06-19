
% This script simulates six sinewaves in two blocks of 3. The signals are dependent within each block.
% It demonstrates MAR modelling to estimate Granger causality between all
% signal pairs

% Parameters
secs=1; % Signal duration
ns=250; % Sampling frequency    
t=[1/ns:1/ns:secs]';    % Time vector
d=6;    % Time series
f1=10;  % Frewuency for first signal group

clear x
% Noise 
dev=1*ones(1,6);

% First block of 3 signals with sinusoid frequency f1 + noise
y=sin(2*pi*f1*t);
y2=sin(2*pi*12.5*t);
x(:,1)=y+dev(1)*randn(size(t));
for i=2:3,
  x(:,i)=y+dev(i)*randn(size(t));
end

% Second block of 3 signals with sinusoid frequency 12.5Hz + noise
for i=4:6,
  x(:,i)=y2+dev(i)*randn(size(t));
end

% Normalise
for i=1:6,
    x(:,i)=x(:,i)/std(x(:,i));  % Variance to 1
    x(:,i)=x(:,i)-mean(x(:,i)); % Zero mean
end

% Fit MAR models
disp('Estimating order of MAR model');
logev=[];
for m=1:5,  % Lags from 1 to 5
    disp(sprintf('Fitting MAR model with %d components',m));
    mar=spm_mar(x,m);
    logev=[logev; mar.fm];  % Store log evidence
end

% Normalise log evidence
logev=logev-min(logev);
figure
bar(logev);
xlabel('Number of time lags');
ylabel('Log Evidence');

% Model order with highest evidence
[tmp, p_sel]=max(logev);

% Fit final MAR model
disp(sprintf('Using MAR(%d) model ..',p_sel));
[mar,y,y_pred]=spm_mar(x,p_sel);

% Calculate Granger causality from fitted MAR model
[G,Psig] = spm_granger (mar);

% Known block structure
disp(' ');
disp('True causality matrix');
[ones(3,3),zeros(3,3);zeros(3,3),ones(3,3)]

% Convert significance matrix to belief matrix
disp('Granger probability matrix:');
Peffect=ones(6,6)-Psig;
Peffect
disp('where ijth entry is our belief that time series i Granger causes j');

% Display Granger causality matrix
figure
imagesc(Peffect);
title('Granger probability matrix');
colormap gray
colorbar
disp(' ');
disp('Inferred Granger causality matrix:');
disp('This is Granger Prob matrix thresholded at 0.95');
Peffect>0.95% Threshold for significant binary matrix





