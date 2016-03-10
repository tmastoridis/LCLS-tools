% This function computes TF using the defined model.
% Error between analytical and measured functions is then computed
% Themis Mastorides/Mastoridis 4/15/2009
% C.H.Rivetta /8/25/2009

function er = fit_TF_pi(x, fit_prm)

% Rescale the parameters to their actual values

x = x.*[fit_prm.fittemp.scale];

% Range of frequencies
w = fit_prm.fittemp.w;

% Carrier frequency
wrf = fit_prm.fittemp.wrf;

switch fit_prm.TF.measurement
    
    case{1} % 'Ana/Dig phase alignment'
        
        % Calculate feedback response
        tauD = fit_prm.digital.tau;
        Hfb = feedback_model(x(1:4), w, tauD, fit_prm);
        
        % Calculate phase response
        Hph = phase_model(x(5:6), w);
        
        % Total response
        Hc = Hfb.*Hph;
        
        % Measured data
        H = fit_prm.TF.H;
        
    case{2} %'Notch identification'
        
        % 'y' are the needed parameters to determine the feedback response.
        % We already know the digital tau, digital gain, and the phase
        % between the two loops, but we are fitting for the analog gain
        
        %         fit_prm.fittemp.rffb.AnalogAttn = prm.rffbI.AnalogAttn;
        % fit_prm.fittemp.rffb.DigitalAttn = prm.rffbI.DigitalAttn;
        % fit_prm.fittemp.rffb.Command_digitGain = prm.rffbI.Command_digitGain;
        % fit_prm.fittemp.rffb.Command_digitIrrFilter= prm.rffbI.Command_digitIrrFilter;
        
        Att_D = fit_prm.fittemp.rffb.DigitalAttn;
        IIR_G = fit_prm.fittemp.rffb.Command_digitGain;
        Gd = (fit_prm.digital.gain/x(1))*(2^IIR_G)*(10^(-Att_D/20));    %gain Dig/Analogue
        
        tauD = fit_prm.digital.tau;
        
        y = [fit_prm.analog.tau x(1) Gd ...
            fit_prm.digital.phase];
        
        % Calculate feedback response
        Hfb = feedback_model(y, w, tauD, fit_prm);
        
        % Calculate phase response
        Hph = phase_model(x(2:3), w);
        
        % Calculate notch response
        Hn = notch_model(x(4:6), w);
        
        % Total response
        Hc = Hfb.*Hph.*Hn;
        
        % Measured data
        H = fit_prm.TF.H;
        
    case{3} % 'OL measurement'
        
        % Calculate cavity response using the set value of 1e5 for Q
        % (Q will be fitted in the second iteration)
        x = [x(1:2) fit_prm.fittemp.Qset x(3:end)];
        Hcav = cavity_model(x(1:3),w,wrf);
        
        % Calculate phase response
        Hph = phase_model(x(4:5), w);
        
        % Total response
        Hc = Hcav.*Hph;
        
        % Measured data (scaled by known responses in fit_loop_pi)
        H = fit_prm.fittemp.H;
        
    case{4} % 'klystron bump compensation'
        
        % Calculate cavity response using the set value of 1e5 for Q (Q
        % will be fitted in the second iteration)
        x = [x(1:2) fit_prm.fittemp.Qset x(3:end)];
        Hcav = cavity_model(x(1:3),w,wrf);
        
        % Calculate phase response
        Hph = phase_model(x(7:8), w);
        
        % Calculate klystron response
        Hb = klystron_model(x(4:6), w, wrf);
        
        % Total response
        Hc = Hcav.*Hph.*Hb;
        
        % Measured data (scaled by known responses in fit_loop_pi)
        H = fit_prm.fittemp.H;
        
    case{5} % 'CL measurement'
        
        % Calculate cavity response using the passed value for Q
        x = [x(1:2) fit_prm.fittemp.Qset x(3:end)];
        Hcav = cavity_model(x(1:3),w,wrf);
        
        % Calculate phase response
        Hph = phase_model(x(4:5), w);
        
        % Total Open loop response
        GH = Hcav.*Hph.*fit_prm.fittemp.Hfb.*fit_prm.fittemp.Hn.*fit_prm.fittemp.Hklys;
        
        % Closed loop response
        Hc = GH./(1-GH);
        
        % Measured data
        H = fit_prm.TF.H;
        
    case{6} % 'klystron bump compensation-measuring before the cavity'
        
        % Calculate klystron response
        Hb = klystron_model(x(1:4), w, wrf);
        
        % Calculate phase response
        Hph = phase_model(x(5:6), w);
        
        % Total response
        Hc = Hph.*Hb;
        
        % Measured data
        H = fit_prm.TF.H;
        
    case{8} % 'OTDFB OL'
        
%         fit_prm.wrf = wrf;
        
        Hc = comb_model(x,w,fit_prm.fittemp);        
        
        H = fit_prm.TF.H;
        
    case{9} % 'OTDFB: OL gain, phase measurement'
        % Calculate feedback response
        tauD = fit_prm.digital.tau;
        
        x = [x(1), x(2), fit_prm.OTFdbk.delay_c];
        
        Hfb = feedback_model(x, w, tauD, fit_prm);
        % Total Open loop response
        Hc = Hfb.*fit_prm.fittemp.Hn.*fit_prm.fittemp.Hcav.*fit_prm.fittemp.Hph.*fit_prm.fittemp.Hklys;
        
        
        H = fit_prm.TF.H;
        
    case{11} % 'OTDFB: CL delay, gain, phase measurement'
        % Calculate feedback response
        tauD = fit_prm.digital.tau;
        
        % Calculate feedback response
        
        Hfb = feedback_model(x, w, tauD, fit_prm);
        % Total Open loop response
        GH = Hfb.*fit_prm.fittemp.Hn.*fit_prm.fittemp.Hcav.*fit_prm.fittemp.Hph.*fit_prm.fittemp.Hklys;
        % Closed loop response
        Hc = GH./(1-GH);
        
        H = fit_prm.TF.H;
        
    case{12} % 'RF FB CL, OTDFB OL: gain, phase, delay measurement'
        
        x = [x(1), x(2), x(3)];
        
        Hcomb = comb_model(x,w,fit_prm.fittemp);  
        % Total Open loop response
        GH = fit_prm.fittemp.Hfb.*fit_prm.fittemp.Hn.*fit_prm.fittemp.Hcav.*fit_prm.fittemp.Hph.*fit_prm.fittemp.Hklys;
        
        Hc = Hcomb.*GH./(1-GH);
        
        H = fit_prm.TF.H;
        
end
%figure(8)
%plot(w/2/pi, 20*log10(abs(H)), w/2/pi, 20*log10(abs(Hc))), grid

%figure(8); plot(w/2/pi, fit_prm.fittemp.Weight.*abs((H)-(Hc)).^2);
% Quadratic error
er = sum(fit_prm.fittemp.Weight.*abs((H)-(Hc)).^2);