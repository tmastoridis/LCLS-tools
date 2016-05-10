function Hc = cavity_model_SS(sys_prm, res_freq, qual_factor, test_freq_arr)

% Cavity Parameters

r_Q = 1036;

fo = res_freq+sys_prm.frf;    % i.e. 1.3e9-1e3
QL = qual_factor; % i.e. 4.12e7

RL = 0.5*r_Q*QL;
wo = 2*pi*fo;

w1_2 = pi*fo/QL;
dw = sys_prm.wrf-wo;

% Electromagnetic cavity model in In/Q frame

Ac = [-w1_2 -dw; dw -w1_2];
Bc = RL*w1_2*[1 0 1 0; 0 1 0 1];
[N,p] = size(Bc);  % N states, p inputs
Cc = eye(2,N);
Dc = zeros(2,4);

cav_C = ss(Ac, Bc, Cc, Dc);
% cav_D = c2d(cav_C, Ts);
% [Ac_D, Bc_D, Cc_D, Dc_D] = ssdata(cav_D);

Hc = freqresp(cav_C, 2*pi*test_freq_arr);
Hc = abs(squeeze(Hc(1,1,:)));
