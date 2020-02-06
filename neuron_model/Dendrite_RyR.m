%% Compute time-domain response of a compartment with ER + RYR's

addpath('../Matlab/mole-master/mole_MATLAB') % path of mole library
%addpath('/Users/pshoemaker/Desktop/MOLE/mole_MATLAB');

dt = 5E-5;   % time step for temporal integration

%% Rate & related constants

% Ca diffusion parameters
DCa = 240; % Ca diffusion coefficient (units um^2.s^-1)
fCa = 1; % fraction to reduce DCa due to intracellular crowding

% Calcium 1st-order (Na exchange) pump constants EFFLUX
kp1 = 75; % units uM.s^-1
kpc1 = 1.8; % units uM

% Calcium leakage rate INFLUX
Klk = 1.5;     % units uM.s^-1

% Calcium 2nd-order (SERCA) pump constants EFFLUX
kp2 = 4400; % units  units uM.s^-1
kpc2 = 0.10; % units uM

% RyR channel parameters INFLUX (Breit & Queisser)
kAneg = 28.8;     % units s^-1
kApos = 1500;     % units uM^-4s^-1
kBneg = 385.9;    % units s^-1
kBpos = 1500;     % units uM^-3s^-1
rkB = kBpos/kBneg;   % ratio of rate constants (for computational use)
kCneg = 0.1;      % units s^-1
kCpos = 1.75;     % units s^-1
% RyR channel flow rate parameter
rCa = 6;          % units 

% rates for calcium buffering (with CalB0=40, gives 5/7 buffered @ equil.)
kcbf = 0.7; % units s^-1.uM^-1 (multiplies [Ca] & [CalB])
kcbb = 0.1;  % units s^-1 (multiplies [CaCalB])

CalB0 = 40; % available Calbindin conc. (uM)   (btot in Breit-Queisser)

VolR = 50; % cytosol/ER volume ratio

CB = 1.2; % Cs threshold (uM) to change computational method for RyR states

%% Scale all rate constants by dt for purposes of temporal integration

% Ca diffusion coefficient: also scale by crowding fractions
DCa = dt*fCa*DCa;

% Calcium 1st-order pump rate
kp1 = dt*kp1;

% Calcium leakage rate (constant
Klk = dt*Klk;
jL = Klk; % set leakage molar increment jL = Klk

% SERCA pump clearance rate
kp2 = dt*kp2;

% RyR receptor state rates
kAneg = kAneg*dt;
kApos = kApos*dt;
kBneg = kBneg*dt;
kBpos = kBpos*dt;
kCneg = kCneg*dt;
kCpos = kCpos*dt;
% Calcium influx rate for open RyR channels
rCa = dt*rCa;

% calcium buffering rates
kcbf = dt*kcbf;
kcbb = dt*kcbb;

%% dendrite morphology

ng = 80; % number of compartments

% location markers for Ca wave computations
m1 = round(0.4*ng);
m2 = m1+1;
m3 = m2+1;

dx = 0.25; % center-to-center spacing of grid compartments (um)

dB = 1/12; % scaling factor for Ca fluxes into cell body

%% Constants for dendritic diffusion computations
k = 2; % Order of Laplacian computation

% for straight dendrites:
L2 = lap(k, ng, dx); % matrix to compute Laplacian (straight dendrites)
LN = robinBC(k, ng, dx, 0, 1); % get kth-order Neumann BC constants
BN2 = -LN(1,2:end)/LN(1,1); % matrix to find boundary value for 0 slope
LN2 = LN(end,:);  % matrix to compute slope @ proximal end
clear LN;

%% Set up initial conditions and history arrays

load equilibrium_state_v1.mat; % quiescent equilibrium initial conditions

% set proximal initial condition: [Ca] = quiescent value
Cb = Cs;
Ce0 = Ce;

% load arrays for dendrite with quiescent parameter values
Cs = Cs*ones(ng+2,1);
Ce = Ce*ones(ng+2,1);
C1 = C1*ones(ng,1);
C2 = C2*ones(ng,1);
O1 = O1*ones(ng,1);
O2 = O2*ones(ng,1);

CaCalB = zeros(ng,1); % bound calcium value

% RATES = [];
% STATES = [];
% CaCalBt = [];

%% Loop on time2

% grids for plots
xgrid1 = [0, 0.5*dx:dx:(ng-0.5)*dx, ng*dx]';

% impose initial condition trigger @ distal end to start response
Cs(1,1) = 2;
Cs(2,1) = 2;

Cst = [];
Cet = [];

Csmax = 0; % to keep track of spatial max Cs

tend = 0.1; % end time
for t = 0:dt:tend % loop on time
    
    % set up figure for cytosolic Ca
    figure(3)
    plot(xgrid1,Cs(:,1))
    xlabel('Distance (uM)');
    ylabel('Concentration (uM)');
    set(gca,'fontsize',15);
    xlim([0 20]);
    ylim([0 10]);
    grid on
%     % set up figure for ER Ca
%     figure(4)
%     plot(xgrid1,Ce(:,1))
%     xlabel('Distance (uM)');
%     ylabel('Concentration (uM)');
%     set(gca,'fontsize',15);
%     xlim([0 20]);
%     ylim([0 300]);
%     grid on
    pause(0.001);

    % Compute Ca wave propagation speed
    [M,I] = max(Cs);
    if ( I>=m1 ) % wait until propagating wave has time to develop
        if M > Csmax
            Csmax = M; % keep running tab of maximum Ca
        end
        figure(3)
        hold on
        plot(xlim, [Csmax,Csmax]); % plot a line at peak value of wave
        hold off
    end
    if ( I==m2 && IP==m1 )
        t0 = t; % mark when wave peak arrives at compartment m2
    end
    if ( I==m3 && IP==m2 && M>0.8*Csmax)
        vel = dx/(t-t0); % compute velocity when peak arrives @ cmprtmt m3
    end
    IP = I; % store location of peak for this iteration
    
    % 'cell body' reservoir Ca
    GCb = LN2*Cs; % grad at proximal end of dendrite
    Cb = Cb - DCa*dB*GCb; % influx (proportional to -grad)
    % efflux: Plasma Membrane(-), SERCA(-)
    Cb = Cb - kp1.*Cb./(kpc1+Cb);
    Cb = Cb - kp2.*Cb.^2./((kpc2+Cb).*Ce0);
 
    % lateral diffusion, cytosolic Ca
    Cs = Cs + DCa*L2*Cs;
    % boundary conditions
    Cs(1,1) = BN2*Cs(2:end,1); % grad(Cs) (i.e. flux) = 0, distal end
    Cs(end,1) = Cb; % boundary condition (proximal end)
    % lateral diffusion, ER Ca
    Ce = Ce + DCa*L2*Ce;
%    % boundary conditions
    Ce(1,1) = BN2*Ce(2:end,1); % grad(Ce) (i.e. flux) = 0, distal end
    Ce(end,1) = BN2*Ce(end-1:-1:1,1); % grad(Ce) = 0, proximal end
    
    % internal Cs for use in computing internal states
    CsM = Cs(2:end-1,1);
    CeM = Ce(2:end-1,1);
    
    % powers of Cs (for Ca-dependent rates)
    Cs2pw = CsM.^2;
    Cs4pw = CsM.^4;
    Cs3pw = CsM.^3;
    
    Open = O1 + O2; % open-channel probability for Ryr
    
    % channel flow increments: Plasma Membrane(-), SERCA(-), RyR(+), total
    jP = kp1.*CsM./(kpc1+CsM);
    jS = kp2.*Cs2pw./((kpc2+CsM).*CeM);
    jR = rCa * Open.* (CeM-CsM);
    jT = jL + jR - jP - jS; % total flux into cytosol
    %RATES = [ RATES; [jP, jL, jS, jR, jT]./dt ];

    % Cytosolic increment
    CsM = CsM + jT;
%    Cst = [ Cst; CsM' ];
    
    % ER increment
    CeM = CeM + VolR*(jS-jR);
%    Cet = [ Cet; CeM' ];
    
    % Ca buffering process
    DCaFB = kcbf.*CsM.*(CalB0-CaCalB);
    DCaBF = kcbb.*CaCalB;
    CsM = CsM + DCaBF - DCaFB;
    CaCalB = CaCalB + DCaFB - DCaBF;
    
    %STATES = [ STATES; C1,O2,C2,O1 ];
    
    % RyR states update -- requires explicit loop due to logical branch
    for k=1:ng
        if CsM(k,1) < CB % full kinetic equations
            C1(k,1) = C1(k,1)+kAneg.*O1(k,1)-kApos.*Cs4pw(k,1).*C1(k,1);
            C2(k,1) = C2(k,1)+kCpos.*O1(k,1)-kCneg.*C2(k,1);
            O2(k,1) = O2(k,1)+kBpos.*Cs3pw(k,1).*O1(k,1)-kBneg.*O2(k,1);
            O1(k,1) = 1 - C1(k,1) - O2(k,1) - C2(k,1);
        else % approximation, for extremely fast Ca-dependent rates
            if C1(k,1) == 0
                dC1O1 = 0;
            else
                dC1O1 = kApos.*Cs4pw(k,1).*C1(k,1) - kAneg.*O1(k,1);
                if dC1O1>C1(k,1)
                    dC1O1 = C1(k,1);
                end
            end
            C1(k,1) = C1(k,1) - dC1O1;
            dC2O1 = kCneg.*C2(k,1) - kCpos.*O1(k,1);
            C2(k,1) = C2(k,1) - dC2O1;
            % following divides increment in proportion to relevant rates
            fO1 = kBneg/(kBpos.*Cs3pw(k,1)+kBneg);
            O1(k,1) = O1(k,1) + fO1*(dC1O1+dC2O1);
            O2(k,1) = O2(k,1) + (1-fO1)*(dC1O1+dC2O1);
        end
    end
    
    % reconstitute arrays w/ boundary values
    Cs = [Cs(1,1);CsM;Cs(end,1)];
    Ce = [Ce(1,1);CeM;Ce(end,1)];
    
end
