clear;clear all;clc;
%v845l.mat
% signal= load('b757l.mat');
% signal= load('b269l.mat');
%% Load Data
signal= load('a311l.mat');

% ECG
ecg_signal = signal.val(1,:);
ecg_signal = ecg_signal - mean(ecg_signal); % Remove DC Offset
raw_ecg = ecg_signal;
% ABP
%% Signal Preprocessing - Remove Noise using Wavelets
% [THR,SORH,KEEPAPP] = ddencmp(IN1,'wv',X)->IN1(Denoise), wv(wavelet),signal
%  wdencmp('gbl',C,L,'wname',N,THR,SORH,KEEPAPP) - performs denoising
[C,L]=wavedec(ecg_signal,8,'db4'); %[approx, detail] as per lecture
[thr,sorh,keepapp]=ddencmp('den','wv',ecg_signal);
cleanecg=wdencmp('gbl',C,L,'db4',8,thr,sorh,keepapp);

% Clean unneeded variables
clear f_y; clear C; clear L;
%% Detrend Data - Remove Baseline Wander
opol = 6; %polynomial degree
ecg_signal = cleanecg;
[p,s,mu] = polyfit((1:numel(ecg_signal)),ecg_signal,opol);
f_y = polyval(p,(1:numel(ecg_signal)),[],mu);
ecg_signal_dt = ecg_signal - f_y;

clear opol; clear s; clear mu;
%% R-Peak Detection in QRS Complex
% 1. Find threshold and locate peaks that are above the threshold
% 2. Look for the pks value in the pre-processed ECG Signal and store the
%    indexes
% 3. Calculate the number of samples between each peak-peak
% 4. Convert this sample number to pk-pk heart rate

maxpk = max(ecg_signal_dt);
R_thresh = maxpk/2;
[pks,locs] = findpeaks(ecg_signal_dt,'MinPeakHeight',R_thresh, 'MinPeakDistance', 5); % Finetune sensitivty (5)
[R_hit, R_hit_idx] = ismember(pks, ecg_signal_dt); % Maps the peaks with the Original ecg_signal
% isequal(R_hit_idx, locs) % Check to see the ismember thing works

for i = 1:numel(R_hit_idx)-1
    pk2pk_samples(i) = R_hit_idx(i+1) - R_hit_idx(i);
end

pk2pk_HR = (pk2pk_samples / 250) * 60;

figure()
plot(pk2pk_HR)
val = (1:numel(ecg_signal_dt));
figure()
plot(val,ecg_signal_dt,val(locs),pks,'or');
xlabel('Samples');
ylabel('Voltage');
title('R-Peak Detection in a True Bradycardia Alarm using the Daubechies Wavelet Transform');


