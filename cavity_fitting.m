function [dw_final, QL_final, Gn_final] = cavity_fitting(sys_prm, noisegain)

close all; figure; hold on;

freq_range = -2e4:0.25:2e4;
w = 2*pi*freq_range;

dw_true = sys_prm.cavity.dw_true;
QL_true = sys_prm.cavity.QL_true;
Gn_true = sys_prm.cavity.Gn_true;

Hc = cavity_model(sys_prm, QL_true, dw_true, Gn_true, w);
Hc_noise=Hc.*(1+noisegain*(randn(size(Hc))+1j*randn(size(Hc))));

QL_guess = 3e8;
dw_guess = -1e4;
Gn_guess = 1.0;

opts = optimset('Display', 'off');

%% FIRST STEP - wo optimization
F = @(x,xdata) cavity_model(sys_prm, QL_guess, x, Gn_guess, xdata);
dw_final = real(lsqcurvefit(F, dw_guess, w, Hc_noise, [], [], opts));

%% SECOND STEP - QL optimization
F = @(x,xdata) cavity_model(sys_prm, x*QL_guess, dw_final, Gn_guess, xdata);
QL_final = QL_guess * real(lsqcurvefit(F, 1, w, Hc_noise, [], [], opts));

%% THIRD STEP - Gain optimization - possibrobably not needed
Gn_final = 1.0;