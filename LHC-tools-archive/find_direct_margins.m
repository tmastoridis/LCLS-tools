% This function computes phase margins of the direct loop model
% estimated by parameters in x. It uses compute_loop(x, prm) to
% compute open-loop transfer function. Returns baseband frequencies
% fm (kHz) and phases pm (deg).
%
% function [fm, pm] = find_direct_margins(x, prm);
%
function [fm, pm] = find_direct_margins(prm)

fstr = inline('abs(dir_gn(x, prm)-1)', 'x', 'prm');

wr = prm.cavity.wr - prm.fittemp.wrf;
wm(1) = fminbnd(fstr, min(prm.fittemp.w), wr, [], prm);
wm(2) = fminbnd(fstr, wr, max(prm.fittemp.w), [], prm);

prm.w=wm;
pm = angle(compute_loop(prm))*180/pi;
fm = wm/2/pi;
