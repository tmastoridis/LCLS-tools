clear all; close all;

sys_prm = model_init();
test_noisegain = 0.0:0.010:0.100;
iters_per_test = 20;

% True output for the cavity_fitting function
x_true = [sys_prm.cavity.dw_true, sys_prm.cavity.QL_true, sys_prm.cavity.Gn_true];

%% Iterative Test

dw_error = zeros(0, 'like', test_noisegain);
QL_error = zeros(0, 'like', test_noisegain);
Gn_error = zeros(0, 'like', test_noisegain);

% Broadcasts to all matlabpool workers that warnings should be disabled
spmd
  warning('off', 'all')
end

for noise_index = 1:length(test_noisegain)
    disp(strcat('Working on noisegain=', num2str(test_noisegain(noise_index))))
    dw_avg = zeros(0, iters_per_test); 
    QL_avg = zeros(0, iters_per_test);
    Gn_avg = zeros(0, iters_per_test);

    % parfor automatically starts the largest possible parpool
    parfor n = 1:iters_per_test
        [dw_iter, QL_iter, Gn_iter] = cavity_fitting(sys_prm, test_noisegain(noise_index));
        dw_avg(n) = dw_iter; 
        QL_avg(n) = QL_iter; 
        Gn_avg(n) = Gn_iter;
    end
    
    % Writes out error values
    error_vals = (x_true - [mean(dw_avg), mean(QL_avg), mean(Gn_avg)]) ./ x_true;
    dw_error(noise_index) = error_vals(1);
    QL_error(noise_index) = error_vals(2);
    Gn_error(noise_index) = error_vals(3);
end