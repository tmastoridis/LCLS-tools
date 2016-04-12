x_true = [-1e3, 1, 6e4];
wrf = 2*pi*4e8;
xdata = 2*pi*(-1e6:1e2:1e6);
noisegain = 1e4;

Hc = cavity_model(x_true, xdata, wrf);
Hc_noise=Hc+noisegain*(randn(size(Hc))+1j*randn(size(Hc)));

x_guess = [-0.5e3, 0.5, 1e4];
%x_guess = [-1, 1, 1];

F = @(x,xdata)cavity_model([x(1), x_guess(2), x_guess(3)], xdata, wrf);
x_fitted = lsqcurvefit(F, x_guess(1), xdata, Hc_noise);
x_guess(1) = x_fitted(1)

F = @(x,xdata)cavity_model([x_guess(1), x(1), x(2)], xdata, wrf);
x_fitted = lsqcurvefit(F, [x_guess(2), x_guess(3)], xdata, Hc_noise);
x_fitted = [x_guess(1), x_fitted(1), x_fitted(2)]

F = @(x,xdata)cavity_model(x, xdata, wrf);

figure
hold on
plot(xdata/2/pi, 20*log10(Hc_noise))
plot(xdata/2/pi, 20*log10(F(x_fitted, xdata)),'-g')
plot(xdata/2/pi, 20*log10(Hc))

figure
plot(xdata/2/pi, F(x_fitted, xdata) ./ Hc)