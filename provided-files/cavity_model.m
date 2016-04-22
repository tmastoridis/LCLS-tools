clear all; close all;

f = -2e3:2e3;

% System Parameters ---------------------

frf = 1.3e9;
wrf = 2*pi*frf;

% % Simulation
% fs = 5e6;
% Ts = 1/fs;

% Cavity Parameters

fo = 1.3e9-1e3;
wo = 2*pi*fo;
QL = 4.12e7;
r_Q = 1036;
len = 1;

w1_2 = pi*fo/QL;

RL = 0.5*r_Q*QL;

% Electromagnetic cavity model in In/Q frame

w1_2 = w1_2;
dw = wrf-wo;

Ac = [-w1_2 -dw; dw -w1_2];
Bc = RL*w1_2*[1 0 1 0; 0 1 0 1];
[N,p] = size(Bc);  % N states, p inputs
Cc = eye(2,N);
Dc = zeros(2,4);

cav_C = ss(Ac, Bc, Cc, Dc);
% cav_D = c2d(cav_C, Ts);
% [Ac_D, Bc_D, Cc_D, Dc_D] = ssdata(cav_D);

H = freqresp(cav_C, 2*pi*f);

figure;
plot(f, abs(squeeze(H(1,1,:))))
