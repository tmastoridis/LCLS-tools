% Notch model
% Themis Mastorides/Mastoridis 4/15/2009

function Hc = notch_model(x,w)

% Passed parameters
wo = x(1); A = x(2); B = x(3); 

% Frequency range (wo, A, B are scaled by 2.pi.4e6, so s is scaled too)
s = j*w/4e6/2/pi;

% Notch response
Hc = (s.^2 + A*s + wo^2)./(s.^2 + B*s + wo^2);