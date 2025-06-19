
% See Brovelli et al. (2004) PNAS 101(26), 9849-9854
% and Geweke (1982) JASA 77 (378), 304-313.
% Note: GEW and PVE to be applied to bivariate data only

close all

% This script uses a MAR model to estimate the Granger causality in the
% frequency domain between two signals

% Noise
noise_dev1=0.1;
noise_dev2=0.01;

% Parameters
secs=1; % Total duration
ns=1000;    % Sample number 
t=[1/ns:1/ns:secs]';
N=length(t);    % Time points
d=2;    % Signal number

% Generate shared oscillatory signal
f1=10;
y=0.5*sin(2*pi*f1.*t)+sin(2*pi*15*t);
y=y/std(y); % Normalise to unit variance

% Delay
delay=50; % Milliseconds
delay_in_samples=ns*delay/1000;

% Signal 1 is a noisy copy of y
y1=y+noise_dev1*randn(N,1);

% Signal 2 is a delayed version of y1
y2=[y1(delay_in_samples:end);zeros(delay_in_samples-1,1)];
y2=y2+noise_dev2*randn(N,1);

y=[y1,y2];  % Combine
disp(sprintf('Signal 2 leads signal 1 by %1.2f ms',delay));

% Plot signals
h=figure;
set(h,'name','Data');
subplot(2,1,1);
plot(t,y1);
title('Signal 1');
subplot(2,1,2);
plot(t,y2);
title('Signal 2');
xlabel('Seconds');

% MAR model
p=10; % Model order
freqs=[0.5:0.5:32]; % Freq range to evaluate
mar = spm_mar(y,p);
mar = spm_mar_spectra (mar,freqs,ns);   % Calculate spectra + Granger measures

% Plot Granger Causality (GEW)
h=figure;
set(h,'name','Granger Causality');
for k=1:d,
    for j=1:d,
        if ~(k==j)
            index=(k-1)*d+j;
            subplot(d,d,index);
            plot(mar.f,mar.gew(:,k,j));
            title(sprintf('From %d to %d',j,k));
        end
    end
end

% Plot Proportion of Variance Explained (PVE)
h=figure;
set(h,'name','Proportion of Variance Explained');
for k=1:d,
    for j=1:d,
        if ~(k==j)
            index=(k-1)*d+j;
            subplot(d,d,index);
            plot(mar.f,mar.pve(:,k,j));
            title(sprintf('From %d to %d',j,k));
            axis([min(mar.f) max(mar.f) 0 1]);
        end
    end
end