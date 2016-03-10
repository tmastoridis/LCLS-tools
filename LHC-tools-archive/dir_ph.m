% dir_ph           - Computes direct-loop phase
%
%   x = dir_ph(w, prm) computes phase of the direct loop model in
%   structure prm at frequency w

function y = dir_ph(w, prm)

prm.fittemp.w = w;
y = angle(compute_loop(prm));
