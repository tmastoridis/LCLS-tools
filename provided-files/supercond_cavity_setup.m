% RF station model LCLS-II 
% Setting the parameters for the cavity model in simulink

clear
close all

%% Parameters

% System Parameters ---------------------

frf = 1.3e9;
wrf = 2*pi*frf;

IbDC = 0.3e-3;

% Cavity Parameters --------------------

fo = 1.3e9;
wo = 2*pi*fo;
QL = 4.12e7;
r_Q = 1036;
len = 1;

w1_2 = pi*fo/QL;

RL = 0.5*r_Q*QL; 

% Mechanical Resonance Parameters ------


fme = [42.72 101.3 172.09 232 281];
Qm = [30.5 848.4 705.3 23.2 108.4];
Klm = [0.105 0.108 0.126 0.092 0.084];

% Solid State Amp (SSA) Parameters

BW = 1e6; 
Pmax = 5e3;

% Controller Paramters -----------------

Kp = 15;
Ki = 3.5e3;


% Delay---------------------------------

Tdelay = 1.7e-6;


% Simulation Parameters ----------------

fs = 5e6;
Ts = 1/fs;
Tfinal = 50e-3;

%% Settings - Operation point

Vin0 = 16.6e6;  
Vq0 = 1;

IbAC = 2*IbDC;

Ii = 1; 
Iq = 0;

% System

% wm = wm_p*sin(2*pi*fm*t);

dw0 = -wrf + wo + 2*pi*(26.88 - 35);


% Electromagnetic cavity model in In/Q frame

w1_2 = w1_2;
dw = 0;

Ac = [-w1_2 -dw; dw -w1_2];
Bc = RL*w1_2*[1 0 1 0; 0 1 0 1];
[N,p] = size(Bc);  % N states, p inputs
Cc = eye(2,N);
Dc = zeros(2,4);

cav_C = ss(Ac, Bc, Cc, Dc);
cav_D = c2d(cav_C, Ts);
[Ac_D, Bc_D, Cc_D, Dc_D] = ssdata(cav_D);

% Mechanical Model ---------------------

wme = 2*pi*fme;
wm_Q = wme./Qm;
LFm = length(fme);

Ami = [0 1; -wme(1)^2 -wme(1)/Qm(1)];
Bmi = [0 1]'*(-Klm(1)*2*pi*wme(1)^2);
%Bm1 = 1e3*[0 1]';
Cmi = [1 0];
Dmi = 0;

Res_LF = ss(Ami, Bmi, Cmi, Dmi);

for i = 2:1:5 %LFm
    
    Ami = [0 1; -wme(i)^2 -wme(i)/Qm(i)];
    Bmi = [0 1]'*(-Klm(i)*2*pi*wme(i)^2);
%Bm1 = 1e3*[0 1]';
    Cmi = [1 0];
    Dmi = 0;
    
    Res_LF = Res_LF + ss(Ami, Bmi, Cmi, Dmi);
end

[Am, Bm, Cm, Dm] = ssdata(Res_LF);

mech_D = c2d(Res_LF, Ts);
[Amt, Bmt, Cmt, Dmt] = ssdata(mech_D);

%Bmt = [Bmt Bmt 0*Bmt];
%Dmt = [Dmt Dmt Dmt];

[Nm,pm] = size(Bm);   % Nm states, pm inputs
[mm,Nm] = size(Cm);   % Nm states, mm outputs

% Complete Cavity model ------


A12 = [-Vq0; Vin0]*Cm;
A21 = 2*Bm*[Vin0 Vq0]/1e12;

Af = [Ac A12; A21 Am];
Bf = [Bc; zeros(Nm,p)];
Cf = [Cc zeros(2,Nm); zeros(mm,N) Cm];
Df = zeros(3,p);

(eig(Af))
eig(Af)/2/pi


% Controller in IN-Q formalism ---------
% PI controller

Api_D = eye(2);
Bpi_D = Ts*eye(2); 
Cpi_D = Ki*eye(2);
Dpi_D = Kp*eye(2);

Ghigh = 0*eye(2);

% Low-pass filters

fl = 0.3e6;   % Hz of noise-limiting low-pass filter
fn = 0.8e6;   % Hz offset of nearby mode to be rejected

sBW = 2*pi*fl;
sBWc = 2*pi*(fl+1i*fn);

zBW = exp(-sBW*Ts);
zBWc = exp(-sBWc*Ts);

Fil1_D = tf(1-zBW, [1 -zBW],Ts);
Fil1_Dc = tf(1-zBW, [1 -zBWc],Ts);

s0 = 2*pi*fn;
z0 = exp(-fn*2*pi*1i*Ts);

Fil1_Dc = (z0 - zBWc)/(z0 - zBW)*Fil1_Dc;

% Implementation en simulink

% real filter

[A1, B1, C1, D1] = ssdata(Fil1_D);

Fil_Dss = ss(A1, B1, C1, D1, Ts);

Fil_IQ_ss = append(Fil_Dss, Fil_Dss);

[A1ss, B1ss, C1ss, D1ss] = ssdata(Fil_IQ_ss);

% complex filter

[A1c, B1c, C1c, D1c] = ssdata(Fil1_Dc);

Fil_DssC = ss(A1c, B1c, C1c, D1c, Ts);

Fil_IQ_ssC = append(Fil_DssC, Fil_DssC);

[A1ssC, B1ssC, C1ssC, D1ssC] = ssdata(Fil_IQ_ssC);

A1ssC = real(A1ssC) + imag(A1ssC)*[0 -1; 1 0];
C1ssC = real(C1ssC) + imag(C1ssC)*[0 -1; 1 0];

Fil_IQ_ssC = ss(A1ssC, B1ssC, C1ssC, D1ssC, Ts);

% Op point - ICs

x = [Vin0;  Vq0];
Ib = [0; 0];

Eacc2 = (x(1)^2 + x(2)^2)/1e12;
xm =  -Am\Bm*Eacc2;

dw = dw0 + Cm*xm

fde = dw/2/pi


% SSA model

Vmax = sqrt(2*50*Pmax);   %Vpeak maximum over 50 ohms

UpperLimit = 10*Vmax;
LowerLimit = 0;

N2 = 700/1e-3/RL;        %Coupling Coeff. estimated now

w_bw = 2*pi*BW;
att_bw = 3;
w_cut = 2*pi*20e6;
att_cut = 15;

[n,wb] = buttord(w_bw, w_cut, att_bw, att_cut, 's');

[B,A] = butter(n,wb,'s');

T11 = tf(B, A);
T11_D = c2d(T11,Ts,'imp');

f = (0:.1:20)*1e6;
[amp,ph] = bode(T11, 2*pi*f);
[amp_D,ph] = bode(T11_D, 2*pi*f);

figure
plot(f,20*log10(squeeze(amp)),f,20*log10(squeeze(amp_D))), grid

Amp_D = [T11 0; 0 T11];

[Assa, Bssa, Cssa, Dssa] = ssdata(T11_D);

Assa = 0; Bssa = 0; Cssa = 0; Dssa = 700;  % Filter bypassed (SSA no freq. resp.)

t = (0:1:100)*Ts;
y_D = lsim(T11_D,ones(size(t)), t);
y = lsim(T11,ones(size(t)), t);

figure
plot(t,y,t,y_D),grid

% Delay (Using Thiran all-pass filter)

delay = thiran(Tdelay,Ts);

[numD, denD] = tfdata(delay,'v');

%% Parameter structure

% Cavity ---------------------------

prm.cavity.a11 = exp(-w1_2*Ts);
prm.cavity.Bssa = Bc_D(:,1:2);
prm.cavity.Bbeam = Bc_D(:,3:4);
prm.cavity.length = len;
prm.cavity.Gs = 1/RL;
prm.cavity.N2 = 1/N2; 

% Mechanical Resonances ------------

prm.mechanical.Amt = Amt;
prm.mechanical.Bmt = Bmt;
prm.mechanical.Cmt = Cmt;
prm.mechanical.Dmt = Dmt;

prm.mechanical.IC = xm;

% SSA ------------------------------

prm.SSA.Assa = Assa;
prm.SSA.Bssa = Bssa;
prm.SSA.Cssa = Cssa;
prm.SSA.Dssa = Dssa;

prm.SSA.Lmax = UpperLimit;
prm.SSA.Lmin = LowerLimit;

% Down Conversion - Attenuation -----

prm.downconv.Attn = 5e-8;
prm.downconv.numDelay = numD;
prm.downconv.denDelay = denD;

% LLRF ------------------------------

prm.LLRF.Api_D = Api_D;
prm.LLRF.Bpi_D = Bpi_D;
prm.LLRF.Cpi_D = Cpi_D;
prm.LLRF.Dpi_D = Dpi_D;

prm.LLRF.Ghigh = Ghigh;

% Filter

prm.LLRF.Afil_D = A1ss;
prm.LLRF.Bfil_D = B1ss;
prm.LLRF.Cfil_D = C1ss;
prm.LLRF.Dfil_D = D1ss;

prm.LLRF.Afil_Dc = A1ssC;
prm.LLRF.Bfil_Dc = B1ssC;
prm.LLRF.Cfil_Dc = C1ssC;
prm.LLRF.Dfil_Dc = D1ssC;