%% Loaded Data
% a142s.mat - True Alarm
% b838s.mat - True Alarm
% b379l.mat - True Alarm [Presentation Data]
% t506s.mat - True Alarm
% v806s.mat - True Alarm [Presentation Data]
% 'b379l.mat'
%% Load Data
clc;clear;clear all;

signal=load('a385l.mat');
ecg_signal=signal.val;
ecg_signal=ecg_signal(1,:);

%% Signal Information
Fs=250; 
T = 1/Fs;
L=length(ecg_signal);
t = (0:L-1)*T;        % Time vector

nyqFreq=Fs/2;
%% Remove DC Component
ecg_mean = mean(ecg_signal);
ecg_noDC = ecg_signal - ecg_mean;
%% Signal Processing - Butterworth Filter
% Baseline wander occurs at < 1 Hz

% Apply Butterworth 8th Order Filter
[b, a] = butter(8, [1 50]/nyqFreq, 'bandpass'); 
ecg_filtered = filter(b,a,ecg_noDC); % Use Signal without DC Offset

%% Signal Processing - Notch Filter (Remove 60 Hz Power Interference)
%d = designfilt('bandstopiir','FilterOrder',2, ...
%               'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
%               'DesignMethod','butter','SampleRate',Fs);
%% FFT Analysis
% Filtered Signal
Y=fft(ecg_filtered);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;

%% Plots
% Plot Raw Signal vs Filtered Signal
figure('Name', 'Raw ECG(Red) vs. Filtered Signal(Blue)')
hold on
title('Raw ECG vs. Filtered ECG')
plot(ecg_signal, 'red')
plot(ecg_filtered, 'blue')
legend('Raw ECG Signal', 'Filtered ECG Signal')
xlabel('Samples')
ylabel('Voltage (V)')
hold off

% Single Sided Spectrum - Filtered Data
plot(f(1:2500),P1(1:2500))
title('Single-Sided Amplitude Spectrum of Filtered ECG Signal')
xlabel('f (Hz)')
ylabel('|ECG(f)|')


%% Extract Heart Rate
[ECG_Max,ECG_i ] = max(P1); % Find max amplitude and corresponding index
bpm_f = f(ECG_i); % Find the matching frequency to the max amplitude
avg_bpm = bpm_f * 60






