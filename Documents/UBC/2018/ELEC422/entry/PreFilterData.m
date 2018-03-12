%% Filter Data

% Filter Attributes
fc=7; % Cut off frequency from SFU article
fn=500; % Nyquist frequency = sample frequency/2;
order = 8; % Filter Order

% Low-pass 8th order butterworth filter
[b a]=butter(order,(fc/fn),'low');

% Filter Data
filtered_data=filter(b,a,data(1,:));

% Plot both filtered and unfiltered data
subplot(2,1,1)
plot(data(1,:))
title('Unfiltered PPG')

subplot(2,1,2)
plot(filtered_data(1,:))
title('Filtered PPG')

%% FFT of Unfiltered and Filtered Data

Fs = 6;            % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = 15000;             % Length of signal
t = (0:L-1)*T;        % Time vector

Y = fft(data(1,:),251);
Pyy = Y.*conj(Y)/251;
f = Fs/251*(0:127);
figure()
plot(f,Pyy(1:128))
title('Power spectral density')
xlabel('Frequency (Hz)')


