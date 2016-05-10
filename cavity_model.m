function Hc = cavity_model(sys_prm, QL, dw, w)

% Cavity Parameters
fo = sys_prm.frf;    % i.e. 1.3e9-1e3
wo = 2*pi*fo;

RL = 0.5*sys_prm.r_Q*QL;

% frequency
s = 1j*w;

Hc = wo*RL/2/QL * (s+wo/2/QL + 1j*dw)./(s.^2 + wo/QL*s +(wo/2/QL)^2 + dw^2); 
% dw = wo - wr