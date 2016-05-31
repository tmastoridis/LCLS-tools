clear all; close all;

f = -2e5:2e5;
s = 2*pi*1j*f;

%% Feedback model

Kp0 = 2.0e3;
x = 0.18;
Kp = x*Kp0;
Ki = 2*pi*2000*x*Kp0;
fl = 0.050e6;   % Hz of noise-limiting low-pass filter
sBW = 2*pi*fl;

Hfb = (Kp + Ki./s)*sBW./(s+sBW);

figure;
subplot(211); plot(f/1e3, 20*log10(abs(Hfb)));
subplot(212); plot(f/1e3, 180/pi*unwrap(angle(Hfb)));

