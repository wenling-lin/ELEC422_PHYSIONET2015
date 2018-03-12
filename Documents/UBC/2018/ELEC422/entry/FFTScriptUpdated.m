%% Constants
%Type of Data
% 1 - PPG
% 2 - AVR
% 3 - ECG
type=3;
%% Filtering Signal
fc=7; % Cut off frequency from SFU article
fn=500; % Nyquist frequency = sample frequency/2;
order = 8; % Filter Order

% Low-pass 8th order butterworth filter
[b a]=butter(order,(fc/fn),'low');

% Filter Data
data = val;
filtered_data=filter(b,a,data(type,:));

%% Plot Unfiltered and Filtered Signal
subplot(2,1,1)
plot(data(type,:))
subplot(2,1,2)
plot(filtered_data)
%% FFT Signal -- Assuming data has loaded

Fs = 250;            % Sampling frequency                    
T = 1/Fs;             % Sampling period     
[row L] = size(filtered_data);
t = (0:L-1)*T;        % Time vector

% Compute the Fourier transform of the signal.
Y_f = fft(filtered_data);

% Unfiltered
Y=fft(data(type,:));

% Compute the two-sided spectrum P2. Then compute the single-sided spectrum P1 based on P2 and the even-valued signal length L
P2_f = abs(Y_f/L);
P1_f = P2_f(1:L/2+1);
P1_f(2:end-1) = 2*P1_f(2:end-1);

% Unfiltered
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

%% Define the frequency domain f and plot the single-sided amplitude spectrum P1.
f = Fs*(0:(L/2))/L;

% Plot Unfiltered
% Zoom in and plot up to certain frequency
figure()
subplot(2,1,1)
plot(f(1:2000),P1(1:2000)) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')


% Plot Filtered
subplot(2,1,2)
plot(f(1:2000),P1_f(1:2000)) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1_f(f)|')


