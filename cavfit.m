clear all; close all; format long;

sys_prm = model_init();
freq_range = -2e3:0.25:2e3;

noisegain = 0.5e9;

fo_true = 1.3e9 - 1e3;
QL_true = 4.12e7;

Hc = cavity_model(sys_prm, fo_true, QL_true, freq_range);
Hc_noise=Hc+noisegain*(randn(size(Hc))+1j*randn(size(Hc)));

fo_guess_base = 1.3e9 - 1e3;
fo_guess_offset = 0;
QL_guess = 4.02e7;

F = @(x, xdata)cavity_model(sys_prm, fo_guess_base + x(1), QL_guess, xdata);
fo_fit_offset = lsqcurvefit(F, fo_guess_offset, freq_range, Hc_noise)
fo_fit = fo_guess_base + fo_fit_offset

F = @(x, xdata)cavity_model(sys_prm, fo_fit, QL_guess, xdata);
QL_fit = lsqcurvefit(F, QL_guess, freq_range, Hc_noise)

figure; hold on;
plot(freq_range, Hc_noise)
plot(freq_range, Hc, '-r')
plot(freq_range, F(fo_fit_offset, freq_range), '-g')
