clearvars;
close all;

[f,Fs] = audioread("original.wav");
T = 1/Fs;
L = length(f);
t = (0:L-1) * T;

N = size(f,1);
figure; stem(t,f);
title('Original: Time-domain'); xlabel('time(seconds)');
%% spectrum
df = Fs /N;
w = (-(N/2):(N/2)-1)*df;
y = fft(f) / N;
y2 = fftshift(y);
figure; plot(w, abs(y2));
title('Original: Amplitude Spectrum'); xlabel('Frequency(Hz)');
