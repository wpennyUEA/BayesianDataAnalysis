
% This script uses AR models to estimate the power spectrum of a time
% series composed of multiple sinusoids

% Parameters
N=200;  % Time points
fs=100; % Sampling frequency
t=[1:1:N]'/fs;  % Time vector

% Three sinusoids
f(1)=16;
f(2)=8;
f(3)=32;
x=zeros(N,1);
for i=1:3,
    x=x+sin(2*pi*f(i)*t);   % Sum
end
x=x+0.1*randn(N,1); % Add small Gaussian noise

% Plot signal
figure
plot(t,x);
xlabel('Seconds');
title('Three sinusoids');

% Fit AR models of increasing order
for p=1:10,
    disp(sprintf('Now fitting model with p=%d coefficients',p));
    ar=spm_ar (x,p,0);  
    logev(p)=ar.fm; % Store model evidence
end

% Normalise
logev=logev-min(logev); 

% Plot model evidence
figure
bar(logev);
ylabel('Log Evidence');
xlabel('Model order');

% Get spectral estimates from model with highest evidence (best model)
[max_log, max_p]=max(logev);
disp(sprintf('AR-%d model has highest evidence',max_p));
ar=spm_ar (x,max_p,0);  % AR model using optimal order

% Power spectrum from AR model
freq=[1:45];    % Range to evaluate
p=spm_ar_freq(ar,freq,fs);

% Plot
figure
plot(freq,p);
xlabel('Frequency');
ylabel('Power');