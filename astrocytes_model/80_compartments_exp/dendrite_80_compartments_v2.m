%% Compute and store states for two dendrites coupled to a cell body
%% this version includes the reaction InP3 -> AA -> (uptake)

addpath('../Matlab/mole-master/mole_MATLAB') % path of mole library
%addpath('/Users/pshoemaker/Desktop/MOLE/mole_MATLAB');

dt = 5E-4;   % time step for temporal integration

ncells = 1;  % single cell terminating dendrite
nmorphs = 1; % no branched dendrites and cell body so we only have morph 2 (branched dendrite)
nstates = 9; % number of relevant state variables

%% Rate & related constants

% rate of InP3 production induced by glutamatergic input
kGI = 40; % units uM.s^-1

DCa = 240; % Ca diffusion coefficient (units um^2.s^-1)
fCa = 1.0; % fraction to reduce DCa due to intracellular crowding
DInP3 = 300; % InP3 diffusion coefficient (units um^2.s^-1)
%%%% fInP3 set to zero => no InP3 diffusion %%%%
fInP3 = 1.0; % fraction to reduce DInP3 due to intracellular crowding

% for [Ca]-dependent production of InP3
kd1f = 5;
% for conversion of InP3 into AA (and thence back into PIP2)
k1b = 2.5;

% for inhibition of InP3 production by PKC
ki2 = 0.0943; % units uM^-1 (multiplies [PKC*])  (k2 in Kang-Othmer)

% for activation of PKC
k4f = 0.6; % units s^-1.uM^-1 (multiplies [Ca])  (k4 in Kang-Othmer)
% for deactivation of PKC*    (k5 in Kang-Othmer)
k4b = 0.5; % units s^-1 

% Calcium infux rate constant for open InP3R channels
kCa = 550;  % units uM.s^-1 (multiplies effective open-channel fraction)
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

% leakage rate
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
%CalB0 = 0;
%load('IC.mat'); % load initial conditions

%% Scale all rate constants by dt for purposes of temporal integration

% Ca & InP3 diffusion coefficients: also scale by crowding fractions
DCa = dt*fCa*DCa;
DInP3 = dt*fInP3*DInP3;

% input-induced InP3 production rate
kGI = dt*kGI;

% PIP2 - InP3, InP3 - AA rates
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

%% cell morphology
% nm(1) = 1; % number of morph #1 (cell body) per cell
% ng(1) = 1; % number of compartments per morph

% NOTE: ng for dendrites constrained by k (next section), must be >= 2*k+1
nm(1) = 1; % number of morph #2 (straight dendrites) per cell $$$$$
ng(1) = 80; % number of compartments per morph

dx = 1; % center-to-center spacing of grid compartments (um)

dB = 1/12; % scaling factor for Ca fluxes into cell body

%% Constants for dendritic diffusion computations
k = 2; % Order of Laplacian computation

m = 80; % number of grid compartments; 2*k+1 = min to support O(2) accuracy
m0 = round(m/5); % # grid compartments to set initial non-0 [Ca] levels

% location markers for Ca wave computations
m1 = round(0.6*m);
m2 = m1+1;
m3 = m2+1;


% for straight dendrites:
L2 = lap(k, ng(1), dx); % matrix to compute Laplacian (straight dendrites)
LN = robinBC(k, ng(1), dx, 0, 1); % get kth-order Neumann BC constants
BN2 = -LN(1,2:end)/LN(1,1); % matrix to find boundary value for 0 slope
LN2 = LN(end,:);  % matrix to compute slope @ proximal end

clear LN;

%% Define cell arrays to hold state & input data

STATES = cell(ncells,nmorphs,nstates);
% index 1: cell number;
% index 2: morphological class:
%   1 => cell body; 2 => straight dendrites; (3 => branched dendrites);
% index 3: states:
%   1 => Ca; 2 => InP3;
%   3 => R; 4 => A; 5 => O; 6 => I2;
%   7 => AA; 8 => CaCalB; 9 => PKCv;

% initial conditions of states for all compartments
IC = [0, 0, 1, 0, 0, 0, 0, 0, 0];
for n = 1:ncells
    for m = 1:nmorphs
        for l = 1:nstates
            % note ng = # compartments, nm = # instances of morph
            if (l==1 || l==2) % diffusive states have ng+2 grid points
                STATES{n,m,l} = IC(l)*ones(ng(m)+2,nm(m));
            else % other states have ng grid points
                STATES{n,m,l} = IC(l)*ones(ng(m),nm(m));
            end
        end
    end
end


% ------------------------------------------------
% This section set's up the input cell
% ------------------------------------------------

INPUTS = cell(ncells,nmorphs);
% index 1: cell number;
% index 2: morphological class:
%   1 => cell body; 2 => straight dendrites; (3 => branched dendrites);

for n = 1:ncells
    for m = 1:nmorphs     
        INPUTS{n,m} = zeros(ng(m),nm(m));
    end
end

% boundary conditions
CaB = 0;
InP3B = 0;

%% Loop on time

% grids for plots
% xgrid1 = [0, dx/2:dx:50-dx/2, 50]';
 xgrid1 = [0, dx/2:dx:80-dx/2, 80]';
% xgrid2 = [dx/2:dx:50-dx/2]';

Camax = 0;

tend = 2.5; % end time ( t(start) = 0)
for t = 0:dt:tend % loop on time
    
    % set up figure
    figure(1)
    plot(xgrid1,STATES{1,1,1}(:,1),'b')
    hold on
    plot(xgrid1,STATES{1,1,2}(:,1),'r')
    grid on
    xlabel('Distance (uM)');
    ylabel('Concentration (uM)');
    set(gca,'fontsize',15);
    xlim([0 82]);
    ylim([0 5]);
    hold off
    
    % ------------------------------------------------
    % This section computes info on the Ca wave.
    % ------------------------------------------------
    
    [M,I] = max(STATES{1,1,1}(:,1));
    if ( I>=m1 ) % wait until propagating wave has time to develop
        if M > Camax
            Camax = M; % keep running tab of maximum Ca
        end
        hold on
        plot(xlim, [Camax,Camax]); % plot a line at peak value of wave
        hold off
    end
    if ( I==m2 && IP==m1 )
        t0 = t; % mark when wave peak arrives at 7*m1
    end
    if ( I==m3 && IP==m2 && M>0.8*Camax)
        vel = dx/(t-t0); % compute velocity when peak arrives at 7*m1-1
    end
    IP = I; % store location of peak for this iteration
    
    
    % ------------------------------------------------
    % This section sets up the Input.
    % ------------------------------------------------
    
    % set up inputs: G-protein pulse in 1st 4 compartments of dendrite 1
    if t<0.1 % 100ms long
        INPUTS{1,1}(1,1) = kGI;
        INPUTS{1,1}(2,1) = kGI;
        INPUTS{1,1}(3,1) = kGI;
        INPUTS{1,1}(4,1) = kGI;
    else
        INPUTS{1,1}(1,1) = 0;
        INPUTS{1,1}(2,1) = 0;
        INPUTS{1,1}(3,1) = 0;
        INPUTS{1,1}(4,1) = 0;
    end    
    
    % ------------------------------------------------
    % This section computes the state updates in individual cells & morphs
    % ------------------------------------------------
    
    for n = 1:ncells % loop on cells
        for m = 1:nmorphs % loop on morphs
            
            %%%% Implement diffusion of Ca and InP3
            
            % retrieve current calcium and InP3 concentrations
            Ca = STATES{n,m,1};
            InP3 = STATES{n,m,2};
            
            % concentration = body concentration, proximal ends
            Ca(end,:) = CaB; % boundary condition (proximal end)
            InP3(end,:) = InP3B; % boundary condition (proximal end)
            % grad (i.e. flux) = 0, distal ends
            Ca(1,:) = BN2*Ca(2:end,:);
            InP3(1,:) = BN2*InP3(2:end,:);
            % lateral diffusion
            Ca = Ca + DCa*L2*Ca;
            InP3 = InP3 + DInP3*L2*InP3;
            % for use in computing internal states
            CaM = Ca(2:end-1,:);
            InP3M = InP3(2:end-1,:);
            
        end
        
        % prevent negative concentrations due to roundoff error
        Ca(Ca<0) = 0;
        InP3(InP3<0) = 0;
                
        %%%% Implement reactions & channel flows
        
        % retrieve states
        R = STATES{n,m,3};
        O = STATES{n,m,4};
        A = STATES{n,m,5};
        I2 = STATES{n,m,6};
        AA = STATES{n,m,7};
        CaCalB = STATES{n,m,8};
        PKCv = STATES{n,m,9};
        
        % retrieve inputs
        Gin = INPUTS{n,m};
        
        % InP3 update
        InP3M = InP3M + (Gin + kd1f*CaM)./(1+ki2*PKCv) - k1b*InP3M;
        
        % AA update
        AA = AA + k1b*InP3M - k1b*AA;
        
        % PKCv update
        PKCv = PKCv + k4f*CaM.*(PKC0-PKCv) - k4b*PKCv;
        
        %   indices into tables of Ca-dependent rates
        iCa = round(5E2*CaM)+1;
        
        % InP3 channel states
        DRO = ph2fA(iCa).*InP3M.*R;
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
        CaM = CaM + kCa*Open;
        
        % Ca influx through putative ARC channels, controlled by AA
        CaM = CaM + kaa*AA;
        
        % Ca efflux due to pumps
        CaMsq = CaM.^2;
        CaM = CaM - kp1.*CaM./(kpc1+CaM) - kp2.*CaMsq./(kpc2sq+CaMsq);
        
        % Ca influx due to leakage
        CaM = CaM + Klk;
        
        % Ca buffering process
        DCaFB = kcbf.*CaM.*(CalB0-CaCalB);
        DCaBF = kcbb.*CaCalB;
        CaM = CaM + DCaBF - DCaFB;
        CaCalB = CaCalB + DCaFB - DCaBF;
        
        % replace arrays in STATES

        STATES{n,m,1} = [Ca(1,:);CaM;Ca(end,:)];
        STATES{n,m,2} = [InP3(1,:);InP3M;InP3(end,:)];

        STATES{n,m,3} = R;
        STATES{n,m,4} = O;
        STATES{n,m,5} = A;
        STATES{n,m,6} = I2;
        STATES{n,m,7} = AA;
        STATES{n,m,8} = CaCalB;
        STATES{n,m,9} = PKCv;
        
    end
    pause(0.0002)
end


%for nn=1:8
%    ICX(nn) = STATES{1,1,nn};
%end
