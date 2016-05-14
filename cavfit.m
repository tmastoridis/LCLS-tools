clear all; close all; format long;

sys_prm = model_init();
freq_range = -2e4:0.25:2e4;
w = 2*pi*freq_range;

noisegain = 0.05;

wo_true = - 1e3 * 2 * pi;
QL_true = 4.12e7;

Hc = cavity_model(sys_prm, QL_true, wo_true, w);
Hc_noise=Hc.*(1+noisegain*(randn(size(Hc))+1j*randn(size(Hc))));

%% FIRST STEP - wo optimization

% Optimization Options (number of iterations etc)
opts = optimset('TolX', 1e-8, 'TolFun', 1e-8);
opts = optimset(opts, 'DiffMaxChange', 1e-2, ...
    'MaxFunEvals', 3e3, 'MaxIter', 1e3);
% Turn 'off' 'iter' at final version (for display purposes only)
opts = optimset(opts,'Display', 'off');

QL_guess = 3e7;
wo_guess = -1e4;

cav_error_wrapper = @(x) cavity_error(sys_prm, [QL_guess, x], w, Hc_noise)
wo_final = fminunc(cav_error_wrapper, wo_guess, opts)
wo_true
wo_error = (wo_final - wo_true) / wo_true

Hc_final = cavity_model(sys_prm, QL_guess, wo_final, w);
figure;
plot(freq_range, 20*log10(Hc_noise), '-r', freq_range, 20*log10(Hc), '-b', freq_range, 20*log10(Hc_final), '-g');

return
%% SECOND STEP - QL optimization

% Optimization Options (number of iterations etc)
opts = optimset('TolX', 1e-8, 'TolFun', 1e-8);
opts = optimset(opts, 'DiffMaxChange', 1e-5, ...
    'MaxFunEvals', 3e3, 'MaxIter', 1e3);
% Turn 'off' 'iter' at final version (for display purposes only)
opts = optimset(opts,'Display', 'off');

cav_error_wrapper = @(x) cavity_error(sys_prm, [x, wo_final], w, Hc_noise)
QL_final = fminunc(cav_error_wrapper, QL_guess, opts)
QL_true
QL_error = (QL_final - QL_true) / QL_true

%% PLOTTING OUTPUT

Hc_final = cavity_model(sys_prm, QL_final, wo_final, w);
figure;
plot(freq_range, 20*lot10(Hc_noise), '-r', freq_range, 20*log10(Hc), '-b', freq_range, 20*log10(Hc_final), '-g');
