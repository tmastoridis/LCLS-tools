function [dw_final, QL_final, Gn_final] = cavity_fitting(sys_prm, noisegain)

freq_range = -2e4:0.25:2e4;
w = 2*pi*freq_range;

dw_true = sys_prm.cavity.dw_true;
QL_true = sys_prm.cavity.QL_true;
Gn_true = sys_prm.cavity.Gn_true;

Hc = cavity_model(sys_prm, QL_true, dw_true, Gn_true, w);
Hc_noise=Hc.*(1+noisegain*(randn(size(Hc))+1j*randn(size(Hc))));

%% FIRST STEP - wo optimization
% Optimization Options (number of iterations etc)
opts = optimset('TolX', 1e-21, 'TolFun', 1e-21);
opts = optimset(opts, 'DiffMaxChange', 1e-8, ...
    'MaxFunEvals', 3e3, 'MaxIter', 2e3);
% Turn 'off' 'iter' at final version (for display purposes only)
opts = optimset(opts,'Display', 'off');

QL_guess = 3e7;
dw_guess = -1e4;
Gn_guess = 1.0;

cav_error_wrapper = @(x) cavity_error(sys_prm, [QL_guess, x, Gn_guess], w, Hc_noise);
dw_final = fminunc(cav_error_wrapper, dw_guess, opts);
dw_error = (dw_final - dw_true) / dw_true;

%% SECOND STEP - QL optimization
opts = optimset('TolX', 1e-17, 'TolFun', 1e-17);
opts = optimset(opts, 'DiffMaxChange', 1e-5, ...
    'MaxFunEvals', 3e3, 'MaxIter', 1e3);
opts = optimset(opts,'Display', 'off');

cav_error_wrapper = @(x) cavity_error(sys_prm, [x, dw_final, Gn_guess], w, Hc_noise);
QL_final = fminunc(cav_error_wrapper, QL_guess, opts);
QL_error = (QL_final - QL_true) / QL_true;

%% THIRD STEP - Gain optimization
opts = optimset('TolX', 1e-15, 'TolFun', 1e-15);
opts = optimset(opts, 'DiffMaxChange', 1e-5, ...
    'MaxFunEvals', 3e3, 'MaxIter', 1e3);
opts = optimset(opts,'Display', 'off');

cav_error_wrapper = @(x) cavity_error(sys_prm, [QL_final, dw_final, x], w, Hc_noise);
Gn_final = fminunc(cav_error_wrapper, Gn_guess, opts);
Gn_error = (Gn_final - Gn_true) / Gn_true;