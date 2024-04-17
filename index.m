%% RLS Adaptive Noise Canceller
%% Gundhada Aman, 510819057
clc;
close all;
clear all;

[rn, fs] = audioread('noisereference.wav'); %g2 reference noise
[x, fs] = audioread('noisyspeech.wav'); %Corrupted Signal
[s,fs] = audioread('denoisedspeech.wav'); %denoised signal from NI

Ts = 1/fs; %Sampling period;
k = 15*fs - 1; %audio recording are 15s long; Max samples.
t = [0:k]*Ts; %Time vector

%Plot Reference Noise Signal
figure, plot(t, rn);
xlabel('Time (t)');
ylabel('Amplitude (V)');
title('Reference Noise');
axis tight
%sound(g2, fs);

%Plot Corrupted speech
figure, plot(t, x),
xlabel('Time (t)');
ylabel('Amplitude (V)');
title('Corrupted Signal');
axis tight
%sound(x, fs);

%Plot Denoised Signal from NI
figure, plot(t, s),
xlabel('Time (t)');
ylabel('Amplitute (V)');
title('Denoised Signal from NI');
axis tight
%sound(s, fs);

%RLS ANC
%RLS FILTER PARAMETERS
M = 7; %Filter length
gamma  = 1; %forgetting factor
gammainv = 1/gamma;
delta = 5;
deltainv = 1/delta;

%Filter Initialization
N = length(x); %number of samples in the input signal signal x
w = zeros(M, 1); %initialize filter Coeefficient
P = deltainv*eye(M); %Inverse correlation matrix
e = zeros(N, 1); % Error Signal

rn = rn(:);
x = x(:);
m = 0; %Number of iterations
r = N - M +1;
%ANC using RLS
for i = M:N
    % Reference input of vector length M
    y = rn(i:-1:i-M+1);
    % Error signal eqn.
    e(i) = x(i) - w'*y;
    %Filter Gain updated
    k = (P*y)/ (gamma + y'*P*y);
    %Inverse correlation matrix updated
    P = (P - k*y'*P)*gammainv;
    %Filter Cooefficients updated using RLS
    w = w + k*e(i);
    fprintf('Progress: %d / %d \n', m , r);
    m = m + 1;
    
    w1(m,:) = w(:, 1);
end
%e = normalize(e, 'range', [-1 1]);

m = min(e);
range = max(e) - m;
e = (e - m)/ range;
e = (e*2) - 1;

figure, 
subplot(3,1,1), plot(s); title('Denoised speech (s) from NI'); axis tight
subplot(3,1,2), plot(x); title('Corrupted Speech (x) '); axis tight,
subplot(3,1,3), plot(e); title('Estimated Speech (e) '); axis tight,
title('ANC using RLS Algorithm');

%Filter Comparisons
figure, 
plot(s, 'r')
hold on
plot(e, '--g')
title('Comparision of Denoised Speech from NI and ANC using RLS Algorithm')
legend('Denoised Speech' , 'ANC output')
axis tight

%Filter Coeffiecients over time
figure,
plot(w1)
title('Adaption of filter coefficient over time');
axis tight























