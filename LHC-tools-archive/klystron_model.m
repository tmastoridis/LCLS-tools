% Klystron Bump Response
% Themis Mastorides/Mastoridis 4/30/2009

function Hc = klystron_model(x,w,wrf) 

% Passed parameters
wrb = x(1); % bump resonance
Gb = x(2); % klystron gain
R = x(3); % bump R
Qr = x(4); % bump Q  

% Frequency range (wrb is scaled by wrf, so is s)
s = j*(w+wrf)/wrf;

% bump response
Hb = (R*wrb/Qr)*s./(s.^2 + (wrb/Qr).*s + wrb^2);
Hc = Gb+Hb;