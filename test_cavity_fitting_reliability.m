clear all; close all; warning('off', 'all');

sys_prm = model_init();

spmd 
  warning('off', 'all')
end

x_true = [sys_prm.cavity.dw_true, sys_prm.cavity.QL_true, sys_prm.cavity.Gn_true];
freq_range = -2e4:0.25:2e4;
ideal = cavity_model(sys_prm, sys_prm.cavity.QL_true, sys_prm.cavity.dw_true, sys_prm.cavity.Gn_true, freq_range);

%% Generate trial data
trial_count = 100;

dw_fit = zeros(0, trial_count);
QL_fit = zeros(0, trial_count);
Gn_fit = zeros(0, trial_count);

parfor n = 1:trial_count
    [dw_fit(n), QL_fit(n), Gn_fit(n)] = cavity_fitting(sys_prm, 0.15);
end

%% Plot trial data
figure; hold on;
plot(freq_range, 20*log10(ideal), '-b')

count_accurate = 0;
count_dwerr = 0;

for n = 1:trial_count
    err = (([dw_fit(n), QL_fit(n), Gn_fit(n)] - x_true) ./ x_true) * 100.0;
    if abs(err(1)) > 1.0
        count_dwerr = count_dwerr + 1;
        plot(freq_range, 20*log10(cavity_model(sys_prm, QL_fit(n), dw_fit(n), Gn_fit(n), freq_range)), '-r')
    else
        count_accurate = count_accurate + 1;
        plot(freq_range, 20*log10(cavity_model(sys_prm, QL_fit(n), dw_fit(n), Gn_fit(n), freq_range)), '-g')
    end
end

[count_accurate, count_dwerr]