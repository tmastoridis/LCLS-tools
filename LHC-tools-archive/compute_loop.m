% compute_loop 
% - This function computes direct loop transfer function
%
% function L = compute_loop(prm)
%
% Computed TF includes
%
% 1. LLRF Digital/Analog loops TF
% 2. Cavity TF
% 3. Klystron TF
% 4. Notch TF
 
function L = compute_loop(fit_prm)
 
wrf = fit_prm.fittemp.wrf;
w = fit_prm.fittemp.w;
 
% 1. LLRF Digital/Analog loops TF-------
        
tauD = fit_prm.digital.tau;
        
Att_A = fit_prm.fittemp.rffb.AnalogAttn;
Ga = fit_prm.analog.gain*(10^(-Att_A/20));                    %gain analog
         
Att_D = fit_prm.fittemp.rffb.DigitalAttn;
IIR_G = fit_prm.fittemp.rffb.Command_digitGain;
Gd = (fit_prm.digital.gain/Ga)*(2^IIR_G)*(10^(-Att_D/20));    %gain Dig/Analogue
 
fb_prm = [fit_prm.analog.tau Ga Gd fit_prm.digital.phase];    %LLRF paramaters
 
Hfb = feedback_model(fb_prm, w, tauD,fit_prm);                        % Trans. Func.
 
% 2. Cavity + phase-Delay TF--------------------------
 
% Cavity paramaters - Trans. Func.
cav_prm = [fit_prm.cavity.wr-wrf fit_prm.cavity.Goo fit_prm.cavity.Q];
Hcav = cavity_model(cav_prm,w,wrf);
 
% Loop phase -Delay paramaters - Trans. Func.
ph_prm = [fit_prm.loop_phase fit_prm.delay.loop];
Hph = phase_model(ph_prm, w);
 
% 3. Klystron TF ------------------------------------
 
klys_prm = [fit_prm.klystron.wr_sim/wrf + 1 fit_prm.klystron.gain_sim ...
            fit_prm.klystron.R_sim fit_prm.klystron.Q_sim]; 
Hklys = klystron_model(klys_prm, w, wrf);
 
% 4. Notch TF
 
notch_prm = [fit_prm.notch.wo/2/pi/4e6 fit_prm.notch.A/2/pi/4e6 ...
     fit_prm.notch.B/2/pi/4e6];
 
Hn = notch_model(notch_prm, w);
 
% Open Loop Transference Function
 
L = Hcav.*Hph.*Hfb.*Hn.*Hklys;
 
 
 

