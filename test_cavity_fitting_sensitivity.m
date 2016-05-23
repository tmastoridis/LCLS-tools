clear all; close all;

N = 500; % number of iterations

sys_prm = model_init();

% True output for the cavity_fitting function
x_true = [sys_prm.cavity.dw_true, sys_prm.cavity.QL_true, sys_prm.cavity.Gn_true];

noisegain = 0.3.*rand(1,N);
dw_guess = 2e1*sys_prm.cavity.dw_true.*(rand(1,N) - 0.5);
QL_guess = 2e2*sys_prm.cavity.QL_true.*rand(1,N);

dw_final = zeros(1, N); QL_final = zeros(1,N); Gn_final = zeros(1,N);

parfor i = 1:N
    [dw_final(i), QL_final(i), Gn_final(i)] = cavity_fitting(sys_prm, noisegain(i), QL_guess(i), dw_guess(i));
end

figure(1)
title('Gain')
plot(noisegain, x_true(3) - Gn_final, 'r.')

figure(2)
title('QL')
plot(QL_guess, x_true(2) - QL_final, 'r.')

figure(3)
title('d\omega')
plot(dw_guess, x_true(1) - dw_final, 'r.')
