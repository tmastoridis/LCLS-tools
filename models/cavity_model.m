% Cavity Response
% Themis Mastorides/Mastoridis 4/15/2009
 
function Hc = cavity_model(x,w,wrf) 
 
RoverQ = 45;
 
% Passed parameters
dw = x(1);      % detuning
Goo = x(2);     % Goo
Q = x(3);       % Q
wr = wrf+dw;    % resonance frequency
 
%Frequency range
s = 1j*(w+wrf);
 
% Cavity response
Hc = Goo*(RoverQ*wr)*s./(s.^2 + (wr/Q).*s + wr^2);
