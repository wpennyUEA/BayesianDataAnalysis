
% This script demonstrates time-resolved spectral analysis using wavelet
% transforms. It shows how wavelets reveal changing frequency content over time

% Parameters
N=200;  % Sample number
fs=100; % Sampling rate
t=[1:1:N]/fs;   % Time vector
freqs=[1:45];   % Freq range to analyse

% Single 10Hz sinusoid
x=sin(2*pi*10*t);
p = spm_wavspec (x,freqs,fs,1); % Wavelet time-frequency analysis

% Plot time-domain signal and spectogram
figure
subplot(2,1,1);
plot(t,x);
title('One sinusoid');
subplot(2,1,2);
imagesc(p); % Power at each time-frequency point
xlabel('Time/Samples');
ylabel('Frequency');

% Two sinusoids - 10Hz and 38Hz
x=x+sin(2*pi*38*t);
p = spm_wavspec (x,freqs,fs);  % Wavelet time-frequency analysis

% Plot combined signal and time-frequency power
figure
subplot(2,1,1);
plot(t,x);
title('Two sinusoids');
subplot(2,1,2);
imagesc(p);
xlabel('Time/Samples');
ylabel('Frequency');

% Chirp
load chirp
p = spm_wavspec (x,freqs,fs);  % Wavelet time-frequency analysis

% Plot chirp signal and spectrogram
figure
subplot(2,1,1);
plot(t,x);
title('Chirp');
subplot(2,1,2);
imagesc(p);
xlabel('Time/Samples');
ylabel('Frequency');