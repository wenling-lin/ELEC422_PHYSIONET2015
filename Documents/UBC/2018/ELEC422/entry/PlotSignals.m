% PPG Signal
subplot(4,1,1)
plot(val(1,:))
title('PPG')

%ABP
subplot(4,1,2)
plot(val(2,:))
title('ABP')

%ECG
subplot(4,1,3)
plot(val(3,:))
title('ECG')

%ECGAVR
subplot(4,1,4)
plot(val(4,:))
title('ECGAVR')