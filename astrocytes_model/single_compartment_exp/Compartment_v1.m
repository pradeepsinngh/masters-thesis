%% Compute and store states for a dendrite coupled to a cell body

addpath('../Matlab/mole-master/mole_MATLAB') % path of mole library
%addpath('/Users/pshoemaker/Desktop/MOLE/mole_MATLAB');

dt = 5E-4;   % time step for temporal integration

%% Rate & related constants

% rate of InP3 production induced by glutamatergic input
kGI = 80; % units uM.s^-1

% for [Ca]-dependent production of InP3
%kd1f = 20; % units s^-1.uM^-1 (multiplies [Ca] & [PIP2])
kd1f = 5;
% for conversion of InP3 back into PIP2
%k1b = 10; % units s^-1 (multiplies [InP3])
k1b = 2.5;

% for inhibition of InP3 production by PKC
ki2 = 0.0943; % units uM^-1 (multiplies [PKC*])  (k2 in Kang-Othmer)
%ki2 = 5; % units uM^-1 (multiplies [PKC*])  (k2 in Kang-Othmer)


% for activation of PKC
k4f = 0.6; % units s^-1.uM^-1 (multiplies [Ca])  (k4 in Kang-Othmer)
% for deactivation of PKC*   (k5 in Kang-Othmer)
k4b = 0.5; % units s^-1

% Calcium infux rate constant for open InP3R channels
kCa = 300;  % units uM.s^-1 (multiplies effective open-channel fraction)
%kCa = 170;  % units uM.s^-1 (multiplies effective open-channel fraction)

% Calcium influx rate for open ARC channels
kaa = 5;  % units uM.s^-1

% Calcium 1st-order pump constants
kp1 = 75; % units uM.s^-1
kpc1 = 1.8; % units uM
% Calcium 2nd-order pump constants
kp2 = 12; % units  units uM.s^-1
kpc2 = 0.10; % units uM
kpc2sq = kpc2^2;

% Ca leakage rate
Klk = 1; % units uM.s^-1

% rates for calcium buffering (with CalB0=40, gives 5/7 buffered @ equil.)
kcbf = 0.7; % units s^-1.uM^-1 (multiplies [Ca] & [CalB])
kcbb = 10; % units s^-1 (multiplies [CaCalB])

% load rate constants / functions for InP3 receptors
load('SD_rates.mat');

% relative dependence of open probablity on O and A states
fO = 0.1;
fA = 0.9;

PKC0 = 1.0; % available PKC concentration (uM)   (KT in Kang-Othmer)
CalB0 = 40; % available Calbindin conc. (uM)   (btot in Breit-Queisser)
%%%% CalB0 set to zero => no calcium buffering
%CalB0 = 0; % set available calbindin to 0 to prevent buffering

%% Scale all rate constants by dt for purposes of temporal integration

% input-induced InP3 production rate
kGI = dt*kGI;

% PIP2 - InP3 rates
kd1f = dt*kd1f;
k1b = dt*k1b;

% Calcium influx rate for open InP3R channels
kCa = dt*kCa;

% Calcium influx rate for open ARC channels
kaa = dt*kaa;

% PKC rates
k4f = dt*k4f;
k4b = dt*k4b;

% Calcium pump/clearance rate
kp1 = dt*kp1;
kp2 = dt*kp2;

% Calcium leakage rate
Klk = dt*Klk;

% calcium buffering rates
kcbf = dt*kcbf;
kcbb = dt*kcbb;

% InP3 receptor state rates -- only those used in reduced model
r5b = dt*r5b;
ph2fA = dt*ph2fA;
ph2bA = dt*ph2bA;
ph4fA = dt*ph4fA;
ph4bA = dt*ph4bA;
ph5fA = dt*ph5fA;

%% Define cell array to hold state data, and fill with initial conditions

% States:
%   1 => Ca; 2 => InP3;
%   3 => R; 4 => A; 5 => O; 6 => I2;
%   7 => CaCalB; 8 => PKCv;

% initial conditions for all compartments
Ca = 0.0;
InP3 = 0;
R = 1;
A = 0;
O = 0;
I2 = 0;
CaCalB = 0;
PKCv = 0;

Cat = [];
InP3t = [];
Rt = [];
At = [];
Ot = [];
I2t = [];
CaCalBt = [];
PKCvt = [];

%% Loop on time

tend = 5; % end time ( t(start) = 0)
for t = 0:dt:tend % loop on time
    
    % G-protein input pulse
    if t<0.1 % 100ms long
        Gin = kGI;
    else
        Gin = 0;
    end
    
    % Update states
    
    % InP3 update : versions with G-protein input
    %InP3 = InP3 + (Gin + kd1f*Ca) - k1b*InP3;
    InP3 = InP3 + (Gin + kd1f*Ca)./(1+ki2*PKCv) - k1b*InP3;
    
    % PKCv update
    PKCv = PKCv + k4f*Ca.*(PKC0-PKCv) - k4b*PKCv;
    
    %   indices into tables of Ca-dependent rates
    iCa = round(5E2*Ca)+1;
    
    % InP3 channel states
    DRO = ph2fA(iCa).*InP3.*R;
    DOR = ph2bA(iCa).*O;
    DOA = ph4fA(iCa).*O;
    DAO = ph4bA(iCa).*A;
    DAI2 = ph5fA(iCa).*A;
    DI2A = r5b.*I2;
    R = R + DOR - DRO;
    O  =  O + DRO + DAO - (DOR + DOA);
    A  =  A + DOA + DI2A - (DAO + DAI2);
    I2 = I2 + DAI2 - DI2A;
    % probability of open state
    Open = (fO*O + fA*A).^4;
    
    % Ca influx through feedback (InP3) channels
    Ca = Ca + kCa*Open;
    
    % Ca influx through putative ARC channels: InP3 = proxy for AA
    Ca = Ca + kaa*InP3;
    
    % Ca efflux due to pumps
    Casq = Ca.^2;
    Ca = Ca - kp1.*Ca./(kpc1+Ca) - kp2.*Casq./(kpc2sq+Casq);
    
    % Ca influx due to leakage
    Ca = Ca + Klk;
    %%%
    
    % Ca buffering process
    DCaFB = kcbf.*Ca.*(CalB0-CaCalB);
    DCaBF = kcbb.*CaCalB;
    Ca = Ca + DCaBF - DCaFB;
    CaCalB = CaCalB + DCaFB - DCaBF;
    
    % prevent negative concentrations due to roundoff error
    %Ca(Ca<0) = 0;
    %InP3(InP3<0) = 0;
    Cat = [Cat,Ca];
    InP3t = [InP3t,InP3];
    Rt = [Rt,R];
    Ot = [Ot,O];
    At = [At,A];
    I2t = [I2t,I2];
    CaCalBt = [CaCalBt,CaCalB];
    PKCvt = [PKCvt,PKCv];
    
end

t=0:dt:tend; % time vector

figure(10)
plot(t,Cat,'b')
hold on
plot(t,InP3t,'r')
plot(t,CaCalBt,'Color',[0 0.5 0]);
%plot(t,PKCvt);

%figure(20)
%plot(t,Rt,'r')
%hold on
%plot(t,Ot,'b')
%plot(t,At,'m')
%plot(t,I2t,'c')
