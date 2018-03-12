%% Loaded Data
% a142s.mat - True Alarm
% b838s.mat - True Alarm
% b379l.mat - True Alarm [Presentation Data]
% t506s.mat - True Alarm

%% Load Data
clc;clear;clear all;

signal=load('b379l.mat');
ecg_signal=signal.val;
ecg_signal=ecg_signal(2,:);

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

% Apply Butterworth 8th Order Filter
%[b, a] = butter(8, [0.49 0.55], 'stop');           
%ecg_filtered = filter(b,a,ecg_signal);


%plot(ecg_filtered, 'blue')
%plot(ecg_signal, 'red')
%hold off
%% Signal Processing - Notch Filter
% notch_filter = designfilt('bandstopiir','FilterOrder',2, ...
%                'HalfPowerFrequency1',0.5,'HalfPowerFrequency2',3, ...
%                'DesignMethod','butter','SampleRate',Fs);
% notchFiltered_ecg = filtfilt(notch_filter,ecg_signal);
%%
% [b, a] = butter(10, [5 7]./(Fs/2), 'stop');
% notchFiltered_ecg = filter(b,a,ecg_noDC);
%% FFT Analysis
%Y=fft(ecg_signal);
%Y=fft(ecg_noDC);
Y=fft(notchFiltered_ecg);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;

%% Plots
figure('Name','ECG Data');
%ECG Raw Signal
subplot(2,1,1);
hold on
title('Raw Unfiltered ECG Signal')
plot(ecg_signal, 'red')
plot(ecg_noDC, 'blue')
xlabel('Samples')
ylabel('Voltage (mV)')
legend('With DC Offset','No DC Offset')

% Single Sided Spectrum
subplot(2,1,2)
plot(f(1:2500),P1(1:2500))
title('Single-Sided Amplitude Spectrum of ECG Signal')
xlabel('f (Hz)')
ylabel('|ECG(f)|')
hold off

% figure('Name', 'Notch Filtered ECG Signal')
% hold on
% title('Notched ECG Signal at 1.5 Hz')
% plot(ecg_signal, 'red');
% plot(notchFiltered_ecg, 'blue');
% xlabel('Samples')
% ylabel('Voltage (mV)')
% legend('Raw ECG Signal','Notched ECG Signal at 1.5 Hz')
% hold off
%% Extract Heart Rate
[ECG_Max,ECG_i ] = max(P1); % Find max amplitude and corresponding index
bpm_f = f(ECG_i); % Find the matching frequency to the max amplitude
avg_bpm = bpm_f * 60






