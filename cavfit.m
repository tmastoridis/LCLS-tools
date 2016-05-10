clear all; close all; format long;

sys_prm = model_init();
freq_range = -2e4:0.25:2e4;
w = 2*pi*freq_range;

noisegain = 1.5e9;

fo_true = - 1e3;
QL_true = 4.12e7;

Hc = cavity_model(sys_prm, QL_true, 2*pi*fo_true, w);

Hc_noise=Hc+noisegain*(randn(size(Hc))+1j*randn(size(Hc)));

plot(freq_range, Hc, '-r', freq_range, Hc_noise, '-r')

%% Noah-Nic
% fo_guess_base = -2e3;
% fo_guess_offset = 0;
% QL_guess = 4.02e7;
% 
% F = @(x, xdata)cavity_model(sys_prm, fo_guess_base + x(1), QL_guess, xdata);
% fo_fit_offset = lsqcurvefit(F, fo_guess_offset, freq_range, Hc_noise)
% fo_fit = fo_guess_base + fo_fit_offset
% 
% F = @(x, xdata)cavity_model(sys_prm, fo_fit, QL_guess, xdata);
% QL_fit = lsqcurvefit(F, QL_guess, freq_range, Hc_noise)
% 
% figure; hold on;
% plot(freq_range, Hc_noise)
% plot(freq_range, Hc, '-r')
% plot(freq_range, F(fo_fit_offset, freq_range), '-g')

%% New

% Optimization Options (number of iterations etc)
opts = optimset('TolX', 1e-8, 'TolFun', 1e-8);
opts = optimset(opts, 'DiffMaxChange', 1e-5, ...
    'MaxFunEvals', 3e3, 'MaxIter', 1e3);
% Turn 'off' 'iter' at final version (for display purposes only)
opts = optimset(opts,'Display', 'off');

x0 = [3e7 -1e4];

% Use function fit_TF_pi to calculate the error function
f = @(x) errorFunction(x, sys_prm, Hc_noise, w);
xFinal = fminunc(f, x0, opts)
HcFinal = cavity_model(sys_prm, xFinal(1), xFinal(2), w);

figure;
plot(freq_range, Hc_noise, '-b', freq_range, HcFinal, '-r')

