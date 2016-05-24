clear all; close all;

load('montecarlo-10000runs-0to0.3noise-10xdwerr-100xqlerr-data.mat')

N = 1e2;
gain_std_data = zeros(1,N);
dw_std_data = zeros(1,N);
QL_std_data = zeros(1,N);
dG = 0.3/N;

QL_data = x_true(2) - QL_final;
dw_data = x_true(1) - dw_final;

for i = 1:N
    a = dG*(i-1);
    b = dG*i;
    ng_pts = find(noisegain >= a & noisegain <= b);
    dw_pts = zeros(1,length(ng_pts));
    QL_pts = zeros(1,length(ng_pts));
    for j = 1:length(ng_pts)
        dw_pts(j) = dw_data(ng_pts(j));
        QL_pts(j) = QL_data(ng_pts(j));
    end
    dw_std_data(i) = std(dw_pts);
    QL_std_data(i) = std(QL_pts);
    gain_std_data(i) = (b+a)/2;
end

figure(1)
hold on
plot(noisegain, abs(dw_data), 'b.')
plot(gain_std_data, dw_std_data, 'r.')

figure(2)
hold on
plot(noisegain, abs(QL_data), 'g.')
plot(gain_std_data, QL_std_data, 'r.')