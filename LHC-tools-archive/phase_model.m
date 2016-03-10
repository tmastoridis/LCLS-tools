% Phase and delay model
% Themis Mastorides/Mastoridis 4/15/2009

function Hc = phase_model(x,w)

% Passed parameters
phi = x(1); % phase
sysDelay = x(2); % delay

% Phase response
Hc = exp(1j*phi*pi/180).*exp(-1j*w*sysDelay);
