% Feedback model (digital and analogue loops incuding the 1Turn-Delay filter)
% C.H.Rivetta 10/25/2009 rivetta@slac.stanford.edu

function Hc = feedback_model(x,w,tauD, fit_prm)


% %     % Passed parameters
%     
% %   tauD :              % digital tau (calculated based on hardware)
%     tauA = x(1);        % analog tau
%     Ga = x(2);          % Ga: analog gain
%     Gd = x(3);          % Gd: digital gain respect to analog gain
%     delta_phi = x(4);   % phase difference between analog and digital paths
%   
% %    % Feedback response
% 
% Hc = -Ga*(Gd*exp(1j*delta_phi*pi/180)./(1+1j*w*tauD) + 1j*w./(1j*w + 1/tauA)); 

if length(x)==3 % 1-Turn-Delay filter 
      
    % 1 Turn-Delay filter Response
    
    H_1TDfdbk = comb_model(x,w,fit_prm.fittemp);
       
    % Digital / Analogue loop parameters
    
%   tauD :                                    % digital tau (calculated based on hardware)
    tauA = fit_prm.fittemp.fb_prm(1) ;        % analog tau
    Ga = fit_prm.fittemp.fb_prm(2);           % Ga: analog gain
    Gd = fit_prm.fittemp.fb_prm(3);           % Gd: digital gain respect to analog gain
    delta_phi = fit_prm.fittemp.fb_prm(4);    % phase difference between analog and digital paths
    
elseif length(x)==4 % Analogue / Digital loop
  
    % Digital / Analogue loop parameters
    
%   tauD :              % digital tau (calculated based on hardware)
    tauA = x(1);        % analog tau
    Ga = x(2);          % Ga: analog gain
    Gd = x(3);          % Gd: digital gain respect to analog gain
    delta_phi = x(4);   % phase difference between analog and digital paths

    % 1 Turn-Delay filter Response
    
    % Check 1-Turn-Delay board status
    if fit_prm.fittemp.OTFdbk.switch == 1; %1-Turn-Delay board connected 
        
        H_1TDfdbk = zeros(size(w));

        %H_1TDfdbk = comb_model(fit_prm.fittemp.OTFdbk_prm,w,fit_prm.fittemp);
    else
        H_1TDfdbk = zeros(size(w));
    end
    
elseif length(x)==7 % Analogue / Digital loop , 1-Turn-Delay filter 
    
    % Digital / Analogue loop parameters
    
%   tauD :              % digital tau (calculated based on hardware)
    tauA = x(1);        % analog tau
    Ga = x(2);          % Ga: analog gain
    Gd = x(3);          % Gd: digital gain respect to analog gain
    delta_phi = x(4);   % phase difference between analog and digital paths
    
   % 1-Turn-Delay filter grouped in a vecto xc = [Gain, phase, delay]
    xc = x(5:7);
    
    H_1TDfdbk = comb_model(xc,w,fit_prm.fittemp);

else
    'length(x) is not defined'
end

%    % Feedback response

Hc = -Ga*(Gd*exp(1j*delta_phi*pi/180)./(1+1j*w*tauD) + 1j*w.*(1+H_1TDfdbk)./(1j*w + 1/tauA)); 
% Hc = -Ga*(1+H_1TDfdbk); % Flat analog loop response, no digital loop or OTFB
