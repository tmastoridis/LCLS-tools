% set_loop_pi      - Computes optimal values for analog and digital loops
%
% This function computes the optimal values for the analog loop gain and
% phase, digital loop gain, phase and one-turn delay based on the measured
% models
%
% function new = set_loop_pi(prm, dpm, cpm)
%
% prm - current (non-optimal) model
% dpm - desired direct loop phase margin (optional, default = 55 degrees)

function new = set_loop_pi(prm)

debug = 0;

dpm = 60;  % Desired Phase Margin  (45-60 deg)

fit_prm = prm.fit_prm;  % To reduce the prm structure size passing to sub-functions
fit_prm.TF.measurement = prm.system.measurement;

fit_prm.fittemp.wrf = prm.system.wrf;
wrf = fit_prm.fittemp.wrf;

% Measured frequencies
%fit_prm.fittemp.w = 2*pi*fit_prm.TF.freq;
%w = fit_prm.fittemp.w;

% Parameters set in the hardware
fit_prm.fittemp.rffb.AnalogAttn = prm.rffbI.AnalogAttn;
fit_prm.fittemp.rffb.DigitalAttn = prm.rffbI.DigitalAttn;
fit_prm.fittemp.rffb.Command_digitGain = prm.rffbI.Command_digitGain;
fit_prm.fittemp.rffb.Command_digitIrrFilter= prm.rffbI.Command_digitIrrFilter;

fit_prm.fittemp.OTFdbk = prm.OTFdbk;
fit_prm.fittemp.h = prm.system.h;


span = 1e6; 

if  max(fit_prm.TF.freq)< span
    fsam = abs(fit_prm.TF.freq(1) - fit_prm.TF.freq(2));
    fit_prm.fittemp.w = 2*pi*(-span:fsam:span);
else
    fit_prm.fittemp.w = 2*pi*[fit_prm.TF.freq];
end

Hc_OLi = compute_loop(fit_prm);
Hc_phase_OLi = [unwrap(angle(Hc_OLi))]; 

if debug
    figure(100);
    subplot(211); plot(fit_prm.w/2/pi, 20*log10(abs(Hc_OLi))), grid
    title('Initial Open-Loop transfer function')
    xlabel('Frequency [Hz.]'), ylabel('Gain (dB)');
    subplot(212); plot(fit_prm.w/2/pi, 180/pi*Hc_phase_OLi), grid
    xlabel('Frequency [Hz.]'), ylabel('Phase (degrees)');
end

new = fit_prm;

% First figure out correct direct loop phase

old_loop_phase = fit_prm.loop_phase;

% For low-Q cavities 
%   fit_prm.loop_phase = 180;
% 
%   fit_prm.w = fit_prm.cavity.wr - fit_prm.wrf;
%   L = compute_loop(fit_prm)

% What the angle should be
%   angle(L)*180/pi
%   new.w(ind_0)
%   new.loop_phase = -angle(L)*180/pi

% For High-Q cavities 

wr = fit_prm.cavity.wr - fit_prm.fittemp.wrf;

k =find(fit_prm.fittemp.w - wr < -2*pi*5e4);
fre_neg = fit_prm.fittemp.w(k) - wr;
phi_neg = Hc_phase_OLi (k);
Pn = polyfit(fre_neg,phi_neg,1);

k =find(fit_prm.fittemp.w - wr > 2*pi*5e4);
fre_pos = fit_prm.fittemp.w(k) - wr;
phi_pos = Hc_phase_OLi (k);
Pp = polyfit(fre_pos,phi_pos,1);

new.loop_phase = fit_prm.loop_phase - (180 + (Pn(2)+Pp(2))*90/pi);
new.loop_phase = rem(new.loop_phase, 360);

xx = (180 + (Pn(2)+Pp(2))*90/pi)
Hc = compute_loop(new);
Hc_phase = [unwrap(angle(Hc))]; 

% correction of phase plot for w = 0Hz, phase [0 -360).

k =find(fit_prm.fittemp.w - wr < -2*pi*5e4);
fre_neg = fit_prm.fittemp.w(k) - wr;
phi_neg = Hc_phase(k);
Pn = polyfit(fre_neg,phi_neg,1);

k =find(fit_prm.fittemp.w - wr > 2*pi*5e4);
fre_pos = fit_prm.fittemp.w(k) - wr;
phi_pos = Hc_phase(k);
Pp = polyfit(fre_pos,phi_pos,1);

wrap = fix((Pn(2)+Pp(2))*90/pi/360);
Hc_phase = Hc_phase - 2*pi*wrap;

if debug
    figure(101);
    subplot(211); plot(new.w/2/pi, 20*log10(abs(Hc))), grid
    title('Open Loop transfer function - Phase corrected')
    xlabel('Frequency [Hz.]'), ylabel('Gain (dB)');
    subplot(212); plot(new.w/2/pi, 180/pi*Hc_phase), grid
    xlabel('Frequency [Hz.]'), ylabel('Phase (degrees)');
end
% Now find direct loop gain for 'dpm' degrees phase margin

% Start from finding current unity gain points (fm in Hz, pm in degrees)

[fm, pm] = find_direct_margins(new);

funstr = inline('abs(dir_ph(x, prm))', 'x', 'prm');

% Margins around the resonance frequency
fm = fm + (new.fittemp.wrf - new.cavity.wr)/2/pi;
fit_prm = new;

% Shift the phase by dpm and find phase zero crossing near -fm
wmin = -2*pi*4e5; 
fit_prm.loop_phase = new.loop_phase + dpm;
low = min(fm)*2*pi;
wm(1) = low;
while (wm(1)==low)
  low = 3000*low;
  wm(1) = fminbnd(funstr, wmin, new.cavity.wr-new.fittemp.wrf, [], fit_prm);
end

% Shift the phase by -dpm and find phase zero crossing near fm
wmax = 2*pi*4e5;
fit_prm.loop_phase = new.loop_phase - dpm;
high = max(fm)*2*pi; wm(2) = high;
while (wm(2)==high)
  high = 2000*high;
  wm(2) = fminbnd(funstr, new.cavity.wr-new.fittemp.wrf, wmax, [], fit_prm);
end

% Find the largest gain
G = max(dir_gn(wm, new));

% Set loop gain so that G=1
new.analog.gain = new.analog.gain/G;
new.digital.gain = new.digital.gain/G;

Hc_OLf = compute_loop(new);

Hc1_phase = [unwrap(angle(Hc_OLf))];
Hc_phase_OLf = Hc1_phase - 2*pi*wrap;

if debug
    figure(102);
    subplot(211); plot(new.w/2/pi, 20*log10(abs(Hc_OLf))), grid
    title('Open Loop transfer function - Gain adjusted')
    xlabel('Frequency [Hz.]'), ylabel('Gain (dB)');
    subplot(212); plot(new.w/2/pi, 180/pi*Hc_phase_OLf), grid
    xlabel('Frequency [Hz.]'), ylabel('Phase (degrees)');
end

% Display the Open Loop Plot-----------------------------------------------   
% determine the range of frequencies -> add the correct legend on the
% x-axis
if max(new.fittemp.w)/2/pi >= 1e6
    freq_label = 'MHz';
    freq_power = 6;
else
    freq_label = 'kHz';
    freq_power = 3;
end
    
figure; clf
subplot(211);
h=plot(new.fittemp.w/2/pi/10^freq_power, 20*log10(abs(Hc_OLi)), 'g', ...
       new.fittemp.w/2/pi/10^freq_power, 20*log10(abs(Hc_OLf)), 'r'); 
set(h(1), 'linewidth', 2);
title('Open Loop transfer function')
legend('Initial', 'Final');
xlabel(['Frequency (',num2str(freq_label),')']);
ylabel('Gain (dB)');
grid;
subplot(212);
h=plot(new.fittemp.w/2/pi/10^freq_power, 180/pi*Hc_phase_OLi, 'g', ...
       new.fittemp.w/2/pi/10^freq_power, 180/pi*Hc_phase_OLf, 'r'); 
set(h(1), 'linewidth', 2);
legend('Initial', 'Final');
xlabel(['Frequency (',num2str(freq_label),')']);
ylabel('Phase (degrees)');
grid;

% CLosed Loop--------------------------------------------------------------

Hc_CLi = Hc_OLi./(1-Hc_OLi);
Hc_phase_CLi = [unwrap(angle(Hc_CLi))];

Hc_CLf = Hc_OLf./(1-Hc_OLf);
Hc_phase_CLf = [unwrap(angle(Hc_CLf))];

%freqcl = new.w/2/pi;
%save close_loop_data Hcl freqcl

if debug
    figure(103)
    subplot(211); plot(new.fittemp.w/2/pi, 20*log10(abs(Hc_CLf))), grid
    title('Closed-Loop Transfer Function')
    xlabel('Frequency [Hz.]'), ylabel('Gain (dB)');
    subplot(212); plot(new.fittemp.w/2/pi, 180/pi*Hc_phase_CLf), grid
    xlabel('Frequency [Hz.]'), ylabel('Phase (degrees)');
end

% Display the Closed Loop Plot-----------------------------------------------   
% determine the range of frequencies -> add the correct legend on the
% x-axis
if max(new.fittemp.w)/2/pi >= 1e6
    freq_label = 'MHz';
    freq_power = 6;
else
    freq_label = 'kHz';
    freq_power = 3;
end
    
figure; clf
subplot(211);
h=plot(new.fittemp.w/2/pi/10^freq_power, 20*log10(abs(Hc_CLi)), 'g', ...
       new.fittemp.w/2/pi/10^freq_power, 20*log10(abs(Hc_CLf)), 'r'); 
set(h(1), 'linewidth', 2);
title('Closed-Loop Transfer Function')
legend('Initial', 'Final');
xlabel(['Frequency (',num2str(freq_label),')']);
ylabel('Gain (dB)');
grid;
subplot(212);
h=plot(new.fittemp.w/2/pi/10^freq_power, 180/pi*Hc_phase_CLi, 'g', ...
       new.fittemp.w/2/pi/10^freq_power, 180/pi*Hc_phase_CLf, 'r'); 
set(h(1), 'linewidth', 2);
legend('Initial', 'Final');
xlabel(['Frequency (',num2str(freq_label),')']);
ylabel('Phase (degrees)');
grid;

%Nyquist plot around the area "real 0, real +1"----------------------------

X_axis_x = [-1 2];
X_axis_y = [0 0];
Y_axis_x = [0 0];
Y_axis_y = [-3 3];
p_one_x = 1;
p_one_y = 0;
p_GM_x = 0.5;
p_GM_y = 0;
L_PM_x = [0 1.5];
L_PM_yn = -tand(dpm)*L_PM_x;
L_PM_yp = tand(dpm)*L_PM_x;
p_PM_x = cosd(dpm);
p_PM_y = sind(dpm);


figure
plot(real(Hc_OLf),imag(Hc_OLf),'.',p_one_x,p_one_y,'k*',X_axis_x ,X_axis_y,'k', Y_axis_x,Y_axis_y,'k',...
    p_GM_x,p_GM_y,'ro',L_PM_x,L_PM_yp,'r',L_PM_x,L_PM_yn,'r',p_PM_x,p_PM_y,'ro',p_PM_x,-p_PM_y,'ro'),grid
title(['Nyquist plot of OL trans. Func.( GH(\omega) ), for PM = ' num2str(dpm) ' deg'])
axis([-1 2 -2 2])
ylabel('IMAG (GH(\omega))'), xlabel('REAL (GH(\omega))');


% Sensibility test---------------------------------------------------------

S = 1./(1-Hc_OLf);  Sn = max(abs(S(1:fix(length(S)/2)))), Sp = max(abs(S(fix(length(S)/2):end)))

figure
h=plot(new.fittemp.w/2/pi/10^freq_power, 20*log10(abs(S)), 'r'); 
set(h(1), 'linewidth', 2);
title('Sensitivity Function')
xlabel(['Frequency (',num2str(freq_label),')']);
ylabel('Sensitivity (dB)');
grid;

    
% Final Numbers

Gdb = 20*log10(1/G);
disp ([' Extra gain G = ' num2str(1/G) '(' num2str(Gdb) ' dB)'] )

rot_phase = rem(new.loop_phase - old_loop_phase, 360);

% set rot_phase  (-180, 180]
if rot_phase  > 180  
    rot_phase  = rot_phase  - 360;
elseif rot_phase <= -180 
    rot_phase  = rot_phase  + 360;
end

disp ([' Rotate the phase ' num2str(rot_phase) ' deg'])

