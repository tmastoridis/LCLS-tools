function Hc = cavity_model(sys_prm, QL, dw, w)
% dw = wo - wr

% Cavity Parameters
fo = sys_prm.frf;
wo = 2*pi*fo;

RL = 0.5*sys_prm.cavity.r_Q*QL;

% frequency
s = 1j*w;

% Cavity transfer function
Hc = wo*RL/2/QL * (s+wo/2/QL + 1j*dw)./(s.^2 + wo/QL*s +(wo/2/QL)^2 + dw^2); 