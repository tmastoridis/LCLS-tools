function [dw_final, QL_final, Gn_final] = cavity_fitting(sys_prm, noisegain)

freq_range = -2e4:0.25:2e4;
w = 2*pi*freq_range;

dw_true = sys_prm.cavity.dw_true;
QL_true = sys_prm.cavity.QL_true;
Gn_true = sys_prm.cavity.Gn_true;

Hc = cavity_model(sys_prm, QL_true, dw_true, Gn_true, w);
Hc_noise=Hc.*(1+noisegain*(randn(size(Hc))+1j*randn(size(Hc))));

%% FIRST STEP - wo optimization
QL_guess = 3e7;
dw_guess = -1e4;
Gn_guess = 1.0;

F = @(x,xdata) cavity_model(sys_prm, QL_guess, x, Gn_guess, xdata);
dw_final = lsqcurvefit(F, dw_guess, w, Hc_noise);

%% SECOND STEP - QL optimization
F = @(x,xdata) cavity_model(sys_prm, x, dw_final, Gn_guess, xdata);
QL_final = lsqcurvefit(F, QL_guess, w, Hc_noise);

%% THIRD STEP - Gain optimization - possibrobably not needed
Gn_final = 1.0;