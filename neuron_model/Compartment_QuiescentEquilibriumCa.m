%% Compute quiescent equilibrium [Ca] in compartment with cytosol & in ER
%  Note: this can be run iteratively,
%  by commenting out the initial conditions after the first run
%  Saving convention: Cs,Ce,C1,C2,O1,O2 in mat file 'equilibrium_state_v#'

dt = 5E-4;   % time step for temporal integration

%% Rate & related constants, v1

% Calcium 1st-order pump constants EFFLUX
kp1 = 75; % units uM.s^-1
kpc1 = 1.8; % units uM

% leakage rate INFLUX
Klk = 1.5;     % units uM.s^-1

% Calcium 2nd-order (SERCA) pump constants EFFLUX
kp2 = 4400; % units  units uM.s^-1
kpc2 = 0.10; % units uM

% RyR channel parameters INFLUX
kAneg = 28.8;     % units s^-1
kApos = 1500;     % units uM^-4s^-1
kBneg = 385.9;    % units s^-1
kBpos = 1500;     % units uM^-3s^-1
kCneg = 0.1;      % units s^-1
kCpos = 1.75;     % units s^-1
rCa = 6; % RyR channel conductance parameter

VolR = 50;  % volume ratio, cytosol to ER

%% initial conditions
Cs = 0.036650138343173;  % initial cytosolic [Ca]
Ce = 2.752319174814520e+02;  % initial ER [Ca]
% initial consitions on RyR states
C1 = 0.998483481216284;
C2 = 0.001423391721516;
O1 = 9.310924525080650e-05;
O2 = 1.781694948757350e-08;
Open = O1 + O2;

%% Loop on time (relaxation)

% arrays to store histories
Cst = [];
Cet = [];
RATES = [];
STATES = [];

tend = 20; % end time ( t(start) = 0)
for t = 0:dt:tend % loop on time

    Cs2pw = Cs.^2;
    Cs4pw = Cs.^4;
    Cs3pw = Cs.^3;
    
    % channel flow rates: Plasma Membrane(-), Leak(+), SERCA(-), RyR(+)
    jP = kp1.*Cs./(kpc1+Cs);
    jL = Klk;
    jS = kp2.*Cs2pw./((kpc2+Cs).*Ce);
    jR = rCa * Open.* (Ce-Cs);
    jT = jL + jR - jP - jS; % total flux into cytosol
    % RATES: 
    RATES = [ RATES; jP, jL, jS, jR, jT ];
    
    % Cytosolic increment
    Cs = Cs + jT*dt;
    Cst = [ Cst, Cs ];
    Cs = max([Cs;0]);
    % ER increment
    Ce = Ce + VolR*(jS-jR)*dt;
    Cet = [ Cet, Ce ];
    
        % RyR states update
    C1 = C1 + (kAneg.*O1 - kApos.*Cs4pw.*C1)*dt;
    O2 = O2 + (kBpos.*Cs3pw.* O1 - kBneg.* O2)*dt;
    C2 = C2 + (kCpos.*O1 - kCneg.*C2)*dt;
    O1 = 1 - C1 - O2 - C2;
    % probability of open state
    Open = O1 + O2;
    STATES = [ STATES; C1,O2,C2,O1 ];

end

t = 0:dt:tend;
figure(1);
plot(t,Cst)
grid on;
figure(2);
plot(t,Cet)
grid on;

save('equilibrium_state_v1.mat','Cs','Ce','C1','C2','O1','O2')