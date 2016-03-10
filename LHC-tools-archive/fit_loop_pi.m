% This script attempts to fit measured data to an idealized transfer
% function.
% Themis Mastoridis, C.H. Rivetta 10/15/2010 (original 4/15/2010)

function prm = fit_loop_pi(prm)

% Assign variables
fit_prm = prm.fit_prm;  % To reduce the prm structure size passing to sub-functions
fit_prm.TF.measurement = prm.system.measurement;
fit_prm.fittemp.wrf = prm.system.wrf;
wrf = fit_prm.fittemp.wrf;
fit_prm.fittemp.Trev = prm.system.Trev;

fit_prm.fittemp.rffb.AnalogAttn = prm.rffbI.AnalogAttn;
fit_prm.fittemp.rffb.DigitalAttn = prm.rffbI.DigitalAttn;
fit_prm.fittemp.rffb.Command_digitGain = prm.rffbI.Command_digitGain;
fit_prm.fittemp.rffb.Command_digitIrrFilter= prm.rffbI.Command_digitIrrFilter;

% Integral Status
if prm.rffbI.Command_digitalFdbkOn && prm.rffbQ.Command_digitalFdbkOn;
    fit_prm.fittemp.rffb.Command_digitalFdbkOn = true;
elseif ~(prm.rffbI.Command_digitalFdbkOn || prm.rffbQ.Command_digitalFdbkOn)
    fit_prm.fittemp.rffb.Command_digitalFdbkOn = false;
else
    if prm.rffbI.Command_digitalFdbkOn
        disp('..')
        disp(' Digital Loop Channel IN is ON / Channel Q is OFF ')
        disp('..')
    elseif prm.rffbQ.Command_digitalFdbkOn
        disp('..')
        disp(' Digital Loop Channel Q is ON / Channel IN is OFF ')
        disp('..')
    end
end

fit_prm.fittemp.OTFdbk = prm.OTFdbk;
fit_prm.fittemp.h = prm.system.h;

% Measured frequencies
fit_prm.fittemp.w = 2*pi*fit_prm.TF.freq;
w = fit_prm.fittemp.w;

% Unwrap phase around wrf (rather than the minimum frequency)
ind_0 = find(w >0 ); ind_0 = ind_0(1);
Hphase = [flipud(unwrap(angle(flipud(fit_prm.TF.H(1:ind_0-1))))); ...
    unwrap(angle(fit_prm.TF.H(ind_0:end)))];

% Set up the weights. These have been determined empirically, by looking
% at the transfer function error (data-estimate) vs frequency.
switch fit_prm.TF.measurement
    case {1,8,9,11} %{'Ana/Dig phase alignment','CL measurement'}
        fit_prm.fittemp.Weight = ones(size(w)); % No weighting
        
    case {3,4} %{'OL measurement','klystron bump compensation','Comb OL','Comb CL'}
        % Weight less around wrf (dominated by cavity)
        ind = find(abs(fit_prm.TF.H(1:ind_0)) == max(abs(fit_prm.TF.H(1:ind_0))));
        detun_est = w(ind(1));
        fit_prm.fittemp.Weight = ((w-detun_est)/max(w)).^2;
        
        % The next two lines are only if measurement 3 will be done with
        % the integral (currently integral is OFF)
        %         ind = find(w<2*pi*1e5 & w>-2*pi*1e5);
        %         fit_prm.fittemp.Weight(ind) = zeros(size(ind));
        
    case {2,5} %{'Notch identification'}
        fit_prm.fittemp.Weight = 1./abs(fit_prm.TF.H);
    case {6}    %{'New Klystron bump Compensation'}
        fit_prm.fittemp.Weight = zeros(size(w));
        wdw = find(w >-2*pi*2e6 & w < 2*pi*5.5e6);
        fit_prm.fittemp.Weight(wdw) = ones(size(wdw));
end

% Initial Conditions
switch fit_prm.TF.measurement
    case{1} % 'Ana/Dig phase alignment'
        
        Ts = prm.system.Ts;
        IIR_fil = prm.rffbI.Command_digitIrrFilter;
        fit_prm.digital.tau = Ts*(2^IIR_fil);   %tau Digital
        
        Att_A = prm.rffbI.AnalogAttn;
        Ga = 3.5*(10^(-Att_A/20));   %Goo (Analog gain, high freq. gain)
        
        Att_D = prm.rffbI.DigitalAttn;
        IIR_G = prm.rffbI.Command_digitGain;
        
        Gd = (0.45/Ga)*(2^IIR_G)*(10^(-Att_D/20));    %gain Dig/Analogue
        
        x0(1)= 100e-6;             %tauA
        x0(2)= Ga;                 %Ga (Analog gain, high freq. gain)
        x0(3)= Gd;                 %Gd gain Dig/Analogue
        x0(4)= 0;                  %rotation in degrees
        x0(5) = 40;                %phase in degrees
        x0(6)= 150e-9;             %delay
        
        % The scale was determined empirically. It is close to the order of
        % magnitude of each variable, but adjusted depending on the partial
        % derivative of the error function with each of the parameters
        % It improves the gradient for each optimization step
        %fit_prm.fittemp.scale = [1e-7 1 50 50 50 1e-6];
        fit_prm.fittemp.scale = [1e-7 .1 5 50 50 1e-6];
        
    case{2} % 'Notch identification'
        
        Att_A = prm.rffbI.AnalogAttn;
        x0(1) = fit_prm.analog.gain*(10^(-Att_A/20));     %gain analog
        
        x0(2) = fit_prm.loop_phase;                         %phase
        x0(3) = 200e-9;                                     %delay
        
        % Notch model is (s.^2 + A*s + wo^2)./(s.^2 + B*s + wo^2);
        %(wo, a, B param. are divided by 2.pi.4e6 to get wo~=1 for fitting purposes)
        
        x0(4) = 2*pi*4e6/2/pi/4e6;                          %wo
        x0(5) = 8e5/2/pi/4e6;                               %A
        x0(6) = 1.1e7/2/pi/4e6;                             %B
        
        % The scale was determined empirically. It is close to the order of
        % magnitude of each variable, but adjusted depending on the partial
        % derivative of the error function with each of the parameters
        % It improves the gradient for each optimization step
        fit_prm.fittemp.scale = [2.8e-2 50 1e-8 1 1e-2 1e-1];
        
    case{3} % 'OL measurement'
        
        % Feedback and Notch parameters are now known. We just pass them in
        % the fittemp structure to be used in the model
        
        tauD = fit_prm.digital.tau;
        
        %If integral OFF, feedback parameters
        
        Att_A = prm.rffbI.AnalogAttn;
        Ga = fit_prm.analog.gain*(10^(-Att_A/20));     %gain analog
        
        fit_prm.fittemp.fb_prm = [fit_prm.analog.tau Ga 0 fit_prm.digital.phase];
        
        %If Integral ON, feedback parameters
        
        %Att_D = prm.rffbI.DigitalAttn;
        %IIR_G = prm.rffbI.Command_digitGain;
        %Gd = (fit_prm.digital.gain/Goo)*(2^IIR_G)*(10^(-Att_D/20));    %gain Dig/Analogue
        
        %fit_prm.fittemp.fb_prm = [fit_prm.analog.tau Ga Gd fit_prm.digital.phase];
        
        % Notch parameters
        fit_prm.fittemp.notch_prm = [fit_prm.notch.wo/2/pi/4e6 fit_prm.notch.A/2/pi/4e6 ...
            fit_prm.notch.B/2/pi/4e6];
        
        % Klystron parameters
        fit_prm.fittemp.klys_prm = [fit_prm.klystron.wr_sim/wrf + 1 fit_prm.klystron.gain_sim ...
            fit_prm.klystron.R_sim fit_prm.klystron.Q_sim];
        
        % Calculate the KNOWN feedback response
        fit_prm.fittemp.Hfb = feedback_model(fit_prm.fittemp.fb_prm, w, tauD,fit_prm);
        
        % Calculate the KNOWN notch response
        fit_prm.fittemp.Hn = notch_model(fit_prm.fittemp.notch_prm, w);
        
        % Calculate the KNOWN SIMPLIFIED klystron response
        fit_prm.fittemp.Hklys = klystron_model(fit_prm.fittemp.klys_prm, w, wrf);
        
        % Scale the measured data by 1/(the known responses), so we only
        % fit the unknown components.
        fit_prm.fittemp.H = fit_prm.TF.H./(fit_prm.fittemp.Hfb.*fit_prm.fittemp.Hn.*fit_prm.fittemp.Hklys);
        
        % Initial Conditions
        x0(1) = detun_est; %2*pi*1e3;   % detuning
        x0(2) = 5e-5;       % Goo
        x0(3) = 0;        % phase
        x0(4) = 660e-9;     % delay
        
        % The scale was determined empirically. It is close to the order of
        % magnitude of each variable, but adjusted depending on the partial
        % derivative of the error function with each of the parameters
        % It improves the gradient for each optimization step
        fit_prm.fittemp.scale = [1e2 1e-5 1e2 1e-7];
        
    case{4} % 'klystron bump compensation'
        
        % Feedback and Notch parameters are now known. We just pass them in
        % the fittemp structure to be used in the model
        
        % Feedback parameters
        fit_prm.fittemp.fb_prm = [fit_prm.digital.tau fit_prm.analog.tau*fit_prm.digital.gain ...
            fit_prm.digital.gain fit_prm.digital.phase];
        
        % Notch parameters
        fit_prm.fittemp.notch_prm = [fit_prm.notch.wo/2/pi/4e6 fit_prm.notch.A/2/pi/4e6 ...
            fit_prm.notch.B/2/pi/4e6];
        
        % Calculate the KNOWN feedback response
        fit_prm.fittemp.Hfb = feedback_model(fit_prm.fittemp.fb_prm, w);
        
        % Calculate the KNOWN notch response
        fit_prm.fittemp.Hn = notch_model(fit_prm.fittemp.notch_prm, w);
        
        % Scale the measured data by 1/(the known responses), so we only
        % fit the unknown components.
        fit_prm.fittemp.H = fit_prm.TF.H./(fit_prm.fittemp.Hfb.*fit_prm.fittemp.Hn);
        
        x0(1) = w(abs(fit_prm.fittemp.H) == max(abs(fit_prm.fittemp.H))); % initial wrf estimate
        x0(2) = 1e-4;               % R/Q
        x0(3) = (wrf+2*pi*4e6)/wrf; % klystron bump resonance frequency (scaled by wrf to get ~1 for fitting purposes)
        x0(4) = .09;                % klystron gain
        x0(5) = 1100;               % klystron bump Q
        x0(6) = -40;                % system phase
        x0(7) = 660e-9;             % system delay
        
        % The scale was determined empirically. It is close to the order of
        % magnitude of each variable, but adjusted depending on the partial
        % derivative of the error function with each of the parameters
        % It improves the gradient for each optimization step
        fit_prm.fittemp.scale = [1e2 1e-2 1e-3 1e-1 1e3 5e2 1e-7]; % [1e2 1e-2 1 1e-1 1e3 5e2 1e-6];
        
    case{5} % 'CL measurement'
        
        x0(1) = -2*pi*1e3;   % detuning
        x0(2) = 5e-5;       % Goo
        x0(3) = 160;          % system phase
        x0(4) = 660e-9;     % system delay
        
        if prm.system.station_name(1) == 'p'
            x0(2) = 5e-6;
        end
        
        % The scale was determined empirically. It is close to the order of
        % magnitude of each variable, but adjusted depending on the partial
        % derivative of the error function with each of the parameters
        % It improves the gradient for each optimization step
        fit_prm.fittemp.scale = [1e3 1e-3 50 1e-7];
        
        % Feedback, Notch and klystron parameters are now known. We just pass them in
        % the fittemp structure to be used in the model
        
        tauD = fit_prm.digital.tau;
        
        Att_A = prm.rffbI.AnalogAttn;
        Ga = fit_prm.analog.gain*(10^(-Att_A/20));     %gain analog
        
        Att_D = prm.rffbI.DigitalAttn;
        IIR_G = prm.rffbI.Command_digitGain;
        Gd = (fit_prm.digital.gain/Ga)*(2^IIR_G)*(10^(-Att_D/20));    %gain Dig/Analogue
        
        fit_prm.fittemp.fb_prm = [fit_prm.analog.tau Ga Gd fit_prm.digital.phase];
        
        % Calculate the KNOWN feedback response
        fit_prm.fittemp.Hfb = feedback_model(fit_prm.fittemp.fb_prm, w, tauD, fit_prm);
        
        % Calculate the KNOWN notch response
        % Notch parameters
        fit_prm.fittemp.notch_prm = [fit_prm.notch.wo/2/pi/4e6 fit_prm.notch.A/2/pi/4e6 ...
            fit_prm.notch.B/2/pi/4e6];
        
        fit_prm.fittemp.Hn = notch_model(fit_prm.fittemp.notch_prm, w);
        
        % Calculate the KNOWN SIMPLIFIED klystron response
        % Klystron parameters
        fit_prm.fittemp.klys_prm = [fit_prm.klystron.wr_sim/wrf + 1 fit_prm.klystron.gain_sim ...
            fit_prm.klystron.R_sim fit_prm.klystron.Q_sim];
        
        fit_prm.fittemp.Hklys = klystron_model(fit_prm.fittemp.klys_prm, w, wrf);
        
    case{6} % 'klystron bump compensation - measuring before the cavity'
        
        x0(1) = (wrf+2*pi*4e6)/wrf; % klystron bump resonance frequency (scaled by wrf to get ~1 for fitting purposes)
        x0(2) = .3;                % klystron gain G
        x0(3) = 1;                 % klystron bump R
        x0(4) = 500;               % klystron bump Q
        x0(5) = 0;                 % system phase
        x0(6) = 400e-9;            % system delay
        
        % The scale was determined empirically. It is close to the order of
        % magnitude of each variable, but adjusted depending on the partial
        % derivative of the error function with each of the parameters
        % It improves the gradient for each optimization step
        fit_prm.fittemp.scale =  [1e-3 1e-1 1 1e3 1e4 1e-7];
        
    case{8,9,11}
        % Set the sub-systems ----------------------------------------
        
        % Feedback, Klystron, Cavity and Notch parameters are now known.
        % We just pass them in the fittemp structure to be used in the
        % model
        
        % Set RFFB modules-------------
        
        Att_A = prm.rffbI.AnalogAttn;
        Ga = fit_prm.analog.gain*(10^(-Att_A/20));     %gain analog
        
        if fit_prm.fittemp.rffb.Command_digitalFdbkOn == 0
            %Integral OFF, feedback parameters
            
            fit_prm.fittemp.fb_prm = [fit_prm.analog.tau Ga 0 fit_prm.digital.phase];
            
        elseif fit_prm.fittemp.rffb.Command_digitalFdbkOn == 1
            %Integral ON, feedback parameters
            
            Att_D = prm.rffbI.DigitalAttn;
            IIR_G = prm.rffbI.Command_digitGain;
            Gd = (fit_prm.digital.gain/Ga)*(2^IIR_G)*(10^(-Att_D/20));    %gain Dig/Analogue
            
            fit_prm.fittemp.fb_prm = [fit_prm.analog.tau Ga Gd fit_prm.digital.phase];
        end
        
        % Set Klystron, Notch -------------
        
        % Notch parameters
        fit_prm.fittemp.notch_prm = [fit_prm.notch.wo/2/pi/4e6 fit_prm.notch.A/2/pi/4e6 ...
            fit_prm.notch.B/2/pi/4e6];
        
        % Klystron parameters
        fit_prm.fittemp.klys_prm = [fit_prm.klystron.wr_sim/wrf + 1 fit_prm.klystron.gain_sim ...
            fit_prm.klystron.R_sim fit_prm.klystron.Q_sim];
        
        % Calculate the KNOWN notch response
        fit_prm.fittemp.Hn = notch_model(fit_prm.fittemp.notch_prm, w);
        
        % Calculate the KNOWN SIMPLIFIED klystron response
        fit_prm.fittemp.Hklys = klystron_model(fit_prm.fittemp.klys_prm, w, wrf);
        
        % Set the cavity and delay model------------------
        
        % Cavity parameters
        fit_prm.fittemp.cav_prm = [fit_prm.cavity.wr - wrf ...
            fit_prm.cavity.Goo fit_prm.cavity.Q];
        
        % Phase - delay parameters
        fit_prm.fittemp.delay_prm = [fit_prm.loop_phase  fit_prm.delay.loop];
        
        % Calculate the KNOWN cavity response
        fit_prm.fittemp.Hcav = cavity_model(fit_prm.fittemp.cav_prm, w, wrf);
        
        % Calculate the KNOWN phase response
        fit_prm.fittemp.Hph = phase_model(fit_prm.fittemp.delay_prm, w);
        
        switch fit_prm.TF.measurement
            case{8}% 'OTDFB: Delay measurement'
                % Initial Conditions -------------------------------------------
                
                x0(1) = 5;             % Gain
                x0(2) = 0;             % phase
                x0(3) = fit_prm.fittemp.Trev;
                fit_prm.fittemp.scale = [1e3 1e1 1e-1];%[1e1 1 1e-6];
                
            case{9}% 'OTDFB: gain, phase open loop'
                
                % Initial Conditions -------------------------------------------
                
                x0(1) = 1;             % Gain
                x0(2) = 0;             % phase
                
                % The scale was determined empirically. It is close to the order of
                % magnitude of each variable, but adjusted depending on the partial
                % derivative of the error function with each of the parameters
                % It improves the gradient for each optimization step
                fit_prm.fittemp.scale = [1e3 1e4]; %[5e2 1e2]*10;
                
            case{11}% 'OTDFB: delay, gain, phase in closed loop'
                
                % Initial Conditions -------------------------------------------
                
                x0(1) = 1;             % Gain
                x0(2) = 0;             % phase
                x0(3) = fit_prm.fittemp.Trev;         % 1T_delay
                
                % The scale was determined empirically. It is close to the order of
                % magnitude of each variable, but adjusted depending on the partial
                % derivative of the error function with each of the parameters
                % It improves the gradient for each optimization step
                fit_prm.fittemp.scale = [1e3 1e4 1e-1];
                
        end
        
end

% For cases {3,4} the fit determines the R/Q value and a second iteration
% follows. The initial value for Q is set to 60e3, except for {5} for which
% it uses the open loop estimate. For {5} the OL or manually entered value
% is used
switch fit_prm.TF.measurement
    case{3,4} % {'OL measurement','klystron bump compensation'}
        fit_prm.fittemp.Qset = 60e3;
    case{5}
        fit_prm.fittemp.Qset = prm.fit_prm.cavity.Q;
end
x0 = x0./[fit_prm.fittemp.scale];

% Optimization Options (number of iterations etc)
opts = optimset('TolX', 1e-8, 'TolFun', 1e-8);
opts = optimset(opts, 'DiffMaxChange', 1e-5, ...
    'MaxFunEvals', 3e3, 'MaxIter', 1e3);
% Turn 'off' 'iter' at final version (for display purposes only)
opts = optimset(opts,'Display', 'off');

% Use function fit_TF_pi to calculate the error function
x = fminunc('fit_TF_pi', x0, opts, fit_prm);

% Rescale the parameters back
x = x.*fit_prm.fittemp.scale;

% separate klystron simplified model fit (case 6)
if fit_prm.TF.measurement == 6
    fit_prm.fittemp.Weight = zeros(size(w));
    wdw = find(w >-2*pi*2e6 & w < 2*pi*2e6);
    fit_prm.fittemp.Weight(wdw) = ones(size(wdw));
    x1 = fminunc('fit_TF_pi', x0, opts, fit_prm);
    
    % Rescale the parameters back
    x1 = x1.*fit_prm.fittemp.scale;
end

% Q fitting. For {3,4} we do a bounded minimization of the model to
% determine the Q value (R/Q was fitted above).
switch fit_prm.TF.measurement
    case{3,4} % {'OL measurement','klystron bump compensation','CL measurement'}
        x = [x(1:2) fit_prm.fittemp.Qset x(3:end)];
        fit_prm.fittemp.x = x;
        
        % Assume that Q values will be between 10k and 210k
        y = fminbnd('fit_Q_pi', 1e4, 2.1e5, opts, fit_prm);  %Qomin = 10K, Qomax = 210K
        x(3) = y;
end


% Plot and Return
switch fit_prm.TF.measurement
    case{1} % 'Ana/Dig phase alignment'
        
        % Determine feedback and phase response
        tauD = fit_prm.digital.tau;
        Hfb = feedback_model(x(1:4), w, tauD, fit_prm);
        Hph = phase_model(x(5:6), w);
        
        % Total response
        Hc = Hfb.*Hph;
        
        % Return parameters. A check of signs of loop phase and digital phase
        % with the appropriate +/- 180 change in phase is done here.
        
        fit_prm.analog.tau = x(1);                 % analog tau
        
        % analog gain conversion
        Att_A = prm.rffbI.AnalogAttn;
        fit_prm.analog.gain = abs(x(2))/(10^(-Att_A/20));
        
        %        fit_prm.digital.tau                       % digital tau
        
        % digital gain conversion
        Att_D = prm.rffbI.DigitalAttn;
        IIR_G = prm.rffbI.Command_digitGain;
        fit_prm.digital.gain = abs(x(3))*abs(x(2))/(2^IIR_G)/(10^(-Att_D/20));
        
        fit_prm.loop_phase = rem(x(5)+(x(2)<0)*180,360);    % loop phase
        
        % set loop_phase between -180 to 180
        if fit_prm.loop_phase > 180
            fit_prm.loop_phase = fit_prm.loop_phase - 360;
        elseif fit_prm.loop_phase <= -180
            fit_prm.loop_phase = fit_prm.loop_phase + 360;
        end
        
        fit_prm.delay.controller = x(6);            % controller delay
        
        % digital phase (between analog and digital paths)
        
        fit_prm.digital.phase = rem(x(4)+((x(3)<0))*180,360);
        
        % set digital_phase between -180 to 180
        if fit_prm.digital.phase > 180
            fit_prm.digital.phase = fit_prm.digital.phase - 360;
        elseif fit_prm.digital.phase <= -180
            fit_prm.digital.phase = fit_prm.digital.phase + 360;
        end
        
        % define the new value that should be set
        prm.setpoint.phases_vcav = rem(prm.setpoint.phases_vcav - fit_prm.digital.phase, 360);
        
        % Plot
        figure; clf
        polar(angle(Hc),abs(Hc), 'g')
        hold on;
        polar(angle(fit_prm.TF.H),abs(fit_prm.TF.H),'r');
        title(['Polar Plot of Measured and Fit Data - Rotation suggested = ' num2str(prm.setpoint.phases_vcav)]);
        
    case{2} % 'Notch identification'
        % Calculate feedback, phase, and notch response
        
        Att_D = prm.rffbI.DigitalAttn;
        IIR_G = prm.rffbI.Command_digitGain;
        Gd = (fit_prm.digital.gain/abs(x(1)))*(2^IIR_G)*(10^(-Att_D/20));    %gain Dig/Analogue
        
        tauD = fit_prm.digital.tau;
        
        y = [fit_prm.analog.tau x(1) Gd ...
            fit_prm.digital.phase];
        
        % Calculate feedback response
        Hfb = feedback_model(y, w, tauD, fit_prm);
        
        Hph = phase_model(x(2:3), w);
        Hn = notch_model(x(4:6), w);
        
        % Total response
        Hc = Hfb.*Hph.*Hn;
        
        % Return parameters. A check of signs with the appropriate +/- 180
        % change in phase is done here.
        
        Att_A = prm.rffbI.AnalogAttn;
        fit_prm.analog.gain = abs(x(1))/(10^(-Att_A/20)); % analog gain conversion
        
        fit_prm.loop_phase = rem(x(2)+(x(1)<0)*180,360);    % LLRF phase
        
        % set loop_phase between -180 to 180
        if fit_prm.loop_phase > 180
            fit_prm.loop_phase = fit_prm.loop_phase - 360;
        elseif fit_prm.loop_phase <= -180
            fit_prm.loop_phase = fit_prm.loop_phase + 360;
        end
        
        fit_prm.delay.controller = x(3);                    % LLRF delay
        
        notch_label = length(prm.fit_prm.notch.DAT); % counter of different capacitor values
        fit_prm.notch.wo_sw(notch_label) = x(4)*2*pi*4e6;   % wo values (array)
        fit_prm.notch.A_sw(notch_label) = x(5)*2*pi*4e6;    % A values (array)
        fit_prm.notch.B_sw(notch_label) = x(6)*2*pi*4e6;    % B values (array)
        
    case{3} % 'OL measurement'
        
        % Calculate cavity and phase responses with fitted values
        Hcav = cavity_model(x(1:3),w,wrf);
        Hph = phase_model(x(4:5), w);
        
        % Total response (including previously known feedback and notch)
        Hc = Hcav.*Hph.*fit_prm.fittemp.Hfb.*fit_prm.fittemp.Hn.*fit_prm.fittemp.Hklys;
        
        % Return parameters
        fit_prm.cavity.wr = x(1) + wrf;         % cavity resonant frequency
        fit_prm.cavity.detune = x(1)/2/pi;      % cavity detune [Hz]
        fit_prm.cavity.Goo = abs(x(2));         % Goo
        fit_prm.cavity.Q = x(3);                % cavity Q
        fit_prm.loop_phase = rem(x(4)+(x(2)<0)*180,360);        % loop phase
        
        % set loop_phase between -180 to 180
        if fit_prm.loop_phase > 180
            fit_prm.loop_phase = fit_prm.loop_phase - 360;
        elseif fit_prm.loop_phase <= -180
            fit_prm.loop_phase = fit_prm.loop_phase + 360;
        end
        
        fit_prm.delay.loop = x(5);              % loop delay
        %        fit_prm.wrf = wrf;                      % RF frequency
        %        fit_prm.w = w;                          % frequency sweep (array)
        
    case{4} % 'klystron bump compensation'
        
        % Calculate cavity, klystron and phase responses with fitted values
        Hcav = cavity_model(x(1:3),w,wrf);
        Hph = phase_model(x(7:8), w);
        Hb = klystron_model(x(4:6), w, wrf);
        
        % Total response (including previously known feedback and notch)
        Hc = Hcav.*Hph.*Hb.*fit_prm.fittemp.Hfb.*fit_prm.fittemp.Hn;
        
        % Return parameters
        fit_prm.cavity.wr = x(1)+wrf; % cavity resonant frequency
        fit_prm.cavity.Goo = abs(x(2)); % Goo
        fit_prm.cavity.Q = x(3); % cavity Q
        fit_prm.klystron.wr = x(4)*wrf; % klystron bump resonant frequency
        fit_prm.klystron.gain = x(5); % klystron gain
        fit_prm.klystron.Q = x(6); % klystorn bump Q
        fit_prm.loop_phase = rem(x(7)+(x(2)<0)*180,360); % loop phase
        fit_prm.delay.loop = x(8); % loop delay
        
    case{5} % 'CL measurement'
        
        % Calculate cavity, klystron and phase responses with fitted values
        Hcav = cavity_model([x(1:2) prm.fit_prm.cavity.Q],w,wrf);
        Hph = phase_model(x(3:4), w);
        
        % Total open loop response (including previously known feedback and notch)
        GH = Hcav.*Hph.*fit_prm.fittemp.Hfb.*fit_prm.fittemp.Hn.*fit_prm.fittemp.Hklys;
        
        % Closed loop response
        Hc = GH./(1-GH);
        
        % Return parameters
        fit_prm.cavity.wr = x(1)+wrf; % cavity resonant frequency
        fit_prm.cavity.Goo = abs(x(2)); % Goo
        fit_prm.cavity.Qcl = prm.fit_prm.cavity.Q; % cavity Q. KEEP OPEN LOOP Q estimate!!!
        fit_prm.loop_phase = rem(x(3)+(x(2)<0)*180,360); % loop phase
        
        % set loop_phase between -180 to 180
        if fit_prm.loop_phase > 180
            fit_prm.loop_phase = fit_prm.loop_phase - 360;
        elseif fit_prm.loop_phase <= -180
            fit_prm.loop_phase = fit_prm.loop_phase + 360;
        end
        
        fit_prm.delay.loop = x(4); % loop delay
        
        
    case{6} % 'klystron bump compensation'
        
        % Calculate cavity, klystron and phase responses with fitted values
        Hb = klystron_model(x(1:4), w, wrf);
        Hph = phase_model(x(5:6), w);
        
        % Total response
        Hc = Hph.*Hb;
        
        % Return parameters
        fit_prm.klystron.wr = (x(1)-1)*wrf; % klystron bump resonant frequency
        fit_prm.klystron.gain = x(2); % klystron gain
        fit_prm.klystron.R = x(3); % klystron bump R
        fit_prm.klystron.Q = x(4); % klystorn Q
        fit_prm.loop_phase = rem(x(5)+(x(2)<0)*180,360); % loop phase
        fit_prm.delay.loop = x(6); % loop delay
        
        fit_prm.klystron.wr_sim = (x1(1)-1)*wrf; % simplified klystron bump resonant frequency
        fit_prm.klystron.gain_sim = x1(2); % simplified klystron gain
        fit_prm.klystron.R_sim = x1(3); % simplified klystron bump R
        fit_prm.klystron.Q_sim = x1(4); % simplified klystorn Q
        fit_prm.klystron.phase_sim = rem(x1(5)+(x1(2)<0)*180,360); % simplified loop phase
        fit_prm.klystron.delay_sim = x1(6); % simplified loop delay
        
    case{8} % 'OTDFB: OL Delay measurement'
        % Calculate rffb transfer function including the "1-turn Delay Filter"
        
        tauD = fit_prm.digital.tau;
        
        % Total response
        Hfb = feedback_model(x, w, tauD, fit_prm);
        % Total Open loop response
        Hc = Hfb.*fit_prm.fittemp.Hn.*fit_prm.fittemp.Hcav.*fit_prm.fittemp.Hph.*fit_prm.fittemp.Hklys;
                
        fit_prm.OTFdbk.Gc = abs(x(1));      % 1-Turn Delay Gain
        fit_prm.OTFdbk.phic  = rem(x(2)+(x(1)<0)*180,360);  % 1-Turn Delay Delay
        % set loop_phase between -180 to 180
        if fit_prm.OTFdbk.phic > 180
            fit_prm.OTFdbk.phic = fit_prm.OTFdbk.phic - 360;
        elseif fit_prm.loop_phase <= -180
            fit_prm.OTFdbk.phic = fit_prm.OTFdbk.phic + 360;
        end
        fit_prm.OTFdbk.delay_c = x(3);      % OTDFB Delay
        
    case{9} % 'OTDFB: OL gain, phase measurement'
        % Calculate rffb transfer function including the "1-turn Delay Filter"
        
        tauD = fit_prm.digital.tau;
        
        % Total response
        x = [x(1), x(2), fit_prm.OTFdbk.delay_c];
        Hfb = feedback_model(x, w, tauD, fit_prm);
        % Total Open loop response
        Hc = Hfb.*fit_prm.fittemp.Hn.*fit_prm.fittemp.Hcav.*fit_prm.fittemp.Hph.*fit_prm.fittemp.Hklys;
        
        % Return parameters
        
        fit_prm.OTFdbk.Gc = abs(x(1));      % 1-Turn Delay Gain
        fit_prm.OTFdbk.phic  = rem(x(2)+(x(1)<0)*180,360);  % 1-Turn Delay Delay
        % set loop_phase between -180 to 180
        if fit_prm.OTFdbk.phic > 180
            fit_prm.OTFdbk.phic = fit_prm.OTFdbk.phic - 360;
        elseif fit_prm.loop_phase <= -180
            fit_prm.OTFdbk.phic = fit_prm.OTFdbk.phic + 360;
        end
        
        
    case{11} % 'OTDFB: CL delay, gain, phase measurement'
        % Calculate rffb transfer function including the "1-turn Delay Filter"
        
        tauD = fit_prm.digital.tau;
        
        % Calculate feedback response
        Hfb = feedback_model(x, w, tauD, fit_prm);
        % Total Open loop response
        GH = Hfb.*fit_prm.fittemp.Hn.*fit_prm.fittemp.Hcav.*fit_prm.fittemp.Hph.*fit_prm.fittemp.Hklys;
        % Closed loop response
        Hc = GH./(1-GH);
        % Return parameters
        
        fit_prm.OTFdbk.Gc = abs(x(1));      % 1-Turn Delay Gain
        fit_prm.OTFdbk.phic  = rem(x(2)+(x(1)<0)*180,360);  % 1-Turn Delay Delay
        % set loop_phase between -180 to 180
        if fit_prm.OTFdbk.phic > 180
            fit_prm.OTFdbk.phic = fit_prm.OTFdbk.phic - 360;
        elseif fit_prm.loop_phase <= -180
            fit_prm.OTFdbk.phic = fit_prm.OTFdbk.phic + 360;
        end
        fit_prm.OTFdbk.delay_c = x(3);      % OTDFB Delay
end

% Unwrap estimated phase around wrf (rather than around the minimum
% frequency)
Hc_phase = [flipud(unwrap(angle(flipud(Hc(1:ind_0-1))))); ...
    unwrap(angle(Hc(ind_0:end)))];

% Plot measured transfer function together with the transfer function
% estimated from the fitted parameters

% determine the range of frequencies -> add the correct legend on the
% x-axis
if max(w)/2/pi >= 1e6
    freq_label = 'MHz';
    freq_power = 6;
else
    freq_label = 'kHz';
    freq_power = 3;
end

figure; clf
subplot(211);
h=plot(w/2/pi/10^freq_power, 20*log10(abs(Hc)), 'g', ...
    w/2/pi/10^freq_power, 20*log10(abs(fit_prm.TF.H)), 'r');
set(h(1), 'linewidth', 2);
legend('Fit', 'Data');
xlabel(['Frequency (',num2str(freq_label),')']);
ylabel('Gain (dB)');
grid;
subplot(212);
h=plot(w/2/pi/10^freq_power, 180/pi*Hc_phase, 'g', ...
    w/2/pi/10^freq_power, 180/pi*Hphase, 'r');
set(h(1), 'linewidth', 2);
xlabel(['Frequency (',num2str(freq_label),')']);
ylabel('Phase (degrees)');
grid;

% Return updated structure
fit_prm = rmfield(fit_prm, 'fittemp');
prm.fit_prm = fit_prm;

