recordName = 'b285l';

%Get all ECG, blood pressure and photoplethysmogram signals
[~,signal,Fs,siginfo]=rdmat(recordName);
alarmResult=1;
description=squeeze(struct2cell(siginfo));
description=description(4,:);

% Resample signal to 125Hz
Fs=Fs(1);
if Fs~=125
    signal=resample(signal,125,Fs);
    Fs=125;
end

h = signal(:,1);

[~,N]=size(h); % h is a 1xN one-dimensional series

nhaf=fix(N/2);

odvn=1;

if nhaf*2==N;
    odvn=0;
end

f=[0:nhaf -nhaf+1-odvn:-1]/N;

Hft=fft(h);

%Compute all frequency domain Gaussians as one matrix

invfk=[1./f(2:nhaf+1)]';

W=2*pi*repmat(f,nhaf,1).*repmat(invfk,1,N);

G=exp((-W.^2)/2); %Gaussian in freq domain

% End of frequency domain Gaussian computation

% Compute Toeplitz matrix with the shifted fft(h)

HW=toeplitz(Hft(1:nhaf+1)',Hft);

% Exclude the first row, corresponding to zero frequency

HW=[HW(2:nhaf+1,:)];

% Compute Stockwell Transform

ST=ifft(HW.*G,[],2); %Compute voice

%Add the zero freq row

st0=mean(h)*ones(1,N);

ST=[st0;ST];


