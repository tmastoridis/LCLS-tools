
% This function computes TF using the defined model.
% Error between analytical and measured functions is then computed
% Themis Mastorides/Mastoridis 4/15/2009

function er = fit_Q_pi(y, fit_prm)

% Fitted parameters from 1st iteration (fit_TF_pi)
x = fit_prm.fittemp.x;

% Digital time constant
tauD = fit_prm.digital.tau;

% carrier frequency
wrf = fit_prm.fittemp.wrf;

% frequency range
w = fit_prm.fittemp.w;
x(3) = y;

switch fit_prm.TF.measurement
    case{3} %'OL measurement'
        
        % empirically it was determined that for a given R/Q ratio 
        % (careful: not a given R!) a variation in Q only affects
        % the open loop transfer function between +/- 10kHz
        
        % use only data between +/- 10 kHz for this fit
        ind = find(w<(x(1)+2*pi*1e4) & w>(x(1)-2*pi*1e4));
        wtemp = w(ind);

        % Compute the phase, notch, and feedback responses based on the
        % parameters from the 1st iteration        
        Hph = phase_model(x(4:5), wtemp);
        Hfb = feedback_model(fit_prm.fittemp.fb_prm, wtemp, tauD, fit_prm);
        Hn = notch_model(fit_prm.fittemp.notch_prm, wtemp);
        % Klystron response correction with SIMPLIFIED model
        Hklys = klystron_model(fit_prm.fittemp.klys_prm(1:4), wtemp, wrf);
            
        % Compute the cavity response with detuning, R/Q from 1st iteration
        % and Q optimized in this function
        Hcav = cavity_model(x(1:3),wtemp,wrf);
                
        % Total response
        Hc = Hcav.*Hph.*Hfb.*Hn.*Hklys;

        % Measured data
        H = fit_prm.TF.H(ind);

    case{4} %'klystron bump compensation'
        
        % empirically it was determined that for a given R/Q ratio 
        % (careful: not a given R!) a variation in Q only affects
        % the open loop transfer function between +/- 10kHz
        
        % use only data between +/- 10 kHz for this fit
        ind = find(w<(x(1)+2*pi*1e4) & w>(x(1)-2*pi*1e4));
        wtemp = w(ind);

        % Compute the phase, klystron, notch, and feedback responses based
        % on the parameters from the 1st iteration
        Hph = phase_model(x(7:8), wtemp);
        Hb = klystron_model(x(4:6), wtemp, wrf);
        Hfb = feedback_model(fit_prm.fittemp.fb_prm, wtemp);
        Hn = notch_model(fit_prm.fittemp.notch_prm, wtemp);
        
        % Compute the cavity response with detuning, R/Q from 1st iteration
        % and Q optimized in this function      
        Hcav = cavity_model(x(1:3),wtemp,wrf);
                
        % Total response
        Hc = Hcav.*Hph.*Hb.*Hfb.*Hn;

        % Measured data
        H = fit_prm.TF.H(ind);

    case{5} %'CL measurement'
        
        % empirically it was determined that for a given R/Q ratio 
        % (careful: not a given R!) a variation in Q only affects
        % the open loop transfer function between +/- 10kHz
        
        % use only data between +/- 10 kHz for this fit
        ind = find(w<(x(1)+2*pi*2e5) & w>(x(1)-2*pi*2e5));
        wtemp = w(ind);
        
        % Compute the cavity response with detuning, R/Q from 1st iteration
        % and Q optimized in this function
        Hcav = cavity_model(x(1:3),wtemp,wrf);
        
        % Compute the phase, notch, and feedback responses based on the
        % parameters from the 1st iteration
        Hph = phase_model(x(4:5), wtemp);
        Hfb = feedback_model(fit_prm.fittemp.fb_prm, wtemp, tauD, fit_prm);
        Hn = notch_model(fit_prm.fittemp.notch_prm, wtemp);
        Hklys = klystron_model(fit_prm.fittemp.klys_prm(1:4), wtemp, wrf);
        
        % Total open loop response
        GH = Hcav.*Hph.*Hfb.*Hn.*Hklys;        
        % Closed loop response
        Hc = GH./(1-GH);

        % Measured data
        H = fit_prm.TF.H(ind);
end

% Quadratic error
er = sum(abs((H)-(Hc)).^2);