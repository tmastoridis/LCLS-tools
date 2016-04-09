x_true = [-1e3, 1, 6e4];
wrf = 2*pi*4e8;
xdata = 2*pi*(-1e6:1e2:1e6);
noisegain = 1e4;

F = @(x,xdata)cavity_model(x,xdata,wrf);

Hc = cavity_model(x_true, xdata, wrf);
Hc_noise=Hc+noisegain*(randn(size(Hc))+1j*randn(size(Hc)));

%x_guess = [-0.5e3, 0.9, 5e4];
%x_guess = [-1e3, 0.1, 50];
x_guess = [-0.5e3, 0.9, 5e4];
x_fitted = lsqcurvefit(F, x_guess, xdata, Hc_noise);

figure
hold on
plot(w/2/pi, 20*log10(Hc_noise))
plot(w/2/pi, 20*log10(F(x_fitted, xdata)),'-g')
plot(w/2/pi, 20*log10(Hc))

%figure 
%hold on
%plot(w/2/pi, 20*log10(Hc))
%plot(w/2/pi, 20*log10(F(x_fitted, xdata)))

x_guess
x_true
abs(x_fitted)
abs((x_true - x_fitted) ./ x_true)

abs(mean((Hc - F(x_fitted, xdata)).^2))