
% This script simulates time-series data from a 5th-order autoregressive process using known coefficients.
% It compares estimated parameters to true ones

% Time parameters
secs=5; % Signal duration
ns=128; % Sampling rate
t=[1/ns:1/ns:secs]; % Time vector
N=length(t);    % Time points

% Noise parameters
noise_var=1^2;  % Varaince of white noise
noise=sqrt(noise_var)*randn(1,N);   % Generate noise

% True AR coefficients
a_true=[-1.8517,1.3741,0.1421,-0.6852,0.3506];

% Generate AR process
y=filter(1,[1,a_true],noise);
y=y(1:N);   % Crop

%Plot simulated signal
figure
plot(t,y);
xlabel('Seconds');
title('Sample of AR-5 process');

% Calculate signal-to-noise ratio
total_variance=std(y)^2;
signal_variance=total_variance-noise_var;
snr=sqrt(signal_variance)/sqrt(noise_var);
disp(sprintf('SNR=%1.3f',snr));

% Model estimation

disp('For model order p=5');
ar=spm_ar (y,5);    % Estimate AR model of order 5

% Compare true and estimated AR coefficients
disp(' ');
disp('True coefficients');
disp(a_true);
disp(' ');
disp('Estimated coefficients');
disp(ar.a_mean');

% Test models from order p=1 to p=10
for p=1:10,
    disp(sprintf('Now fitting model with p=%d coefficients',p));
    ar=spm_ar (y,p,0);
    logev(p)=ar.fm; % Model evidence
end

% Normalise log evidences
logev=logev-min(logev);

% Plot model evidence
figure
bar(logev);
ylabel('Log Evidence');
xlabel('Model order');