
%% Compute and store states for two dendrites coupled to a cell body
% this version includes the reaction InP3 -> AA -> (uptake)

addpath('../Matlab/mole-master/mole_MATLAB') % path of mole library
%addpath('/Users/pshoemaker/Desktop/MOLE/mole_MATLAB');

dt = 5E-4;   % time step for temporal integration
tend = 5; % end time ( t(start) = 0)
itend = floor(tend/dt)+1;

ncells = 2;  % single cell terminating dendrite
nmorphs = 2; % no branched dendrites, so we only have morph 1 and morp 2
nstates = 9; % number of relevant state variables

%% Rate & related constants

% rate of InP3 production induced by glutamatergic input
kGI = 100; % units uM.s^-1
%kGI = 50; % units uM.s^-1

DCa = 240; % Ca diffusion coefficient (units um^2.s^-1)
% fCa = 1.0; % fraction to reduce DCa due to intracellular crowding
DInP3 = 300; % InP3 diffusion coefficient (units um^2.s^-1)
%%%% fInP3 set to zero => no InP3 diffusion %%%%
% fInP3 = 1.0; % fraction to reduce DInP3 due to intracellular crowding

% for [Ca]-dependent production of InP3
% kd1f = 5;
% for conversion of InP3 into AA (and thence back into PIP2)
k1b = 2.5;

% for inhibition of InP3 production by PKC
% ki2 = 0.0943; % units uM^-1 (multiplies [PKC*])  (k2 in Kang-Othmer)

% for activation of PKC
k4f = 0.6; % units s^-1.uM^-1 (multiplies [Ca])    (k4 in Kang-Othmer)
% for deactivation of PKC*    (k5 in Kang-Othmer)
k4b = 0.5; % units s^-1 

% Calcium infux rate constant for open InP3R channels
%kCa = 570;  % units uM.s^-1 (multiplies effective open-channel fraction)

% Calcium influx rate for open ARC channels
% kaa = 4;  % units uM.s^-1

% Calcium 1st-order pump constants
kp1 = 75; % units uM.s^-1
kpc1 = 1.8; % units uM
% Calcium 2nd-order pump constants
kp2 = 12; % units  units uM.s^-1
kpc2 = 0.10; % units uM
kpc2sq = kpc2^2;

% rate constant at interconnections
% rica = 500;
% rinp3 = 500;

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
%load('IC.mat'); % load initial conditions

%% Exp with these parameters:

fCa = 0.3;
fInP3 = 0.7;
kCa = 600;
kaa = 6;
kd1f = 6;
ki2 = 0.0943;
rinp3 = 500;


rica =500;

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

rica = dt*rica;
rinp3 = dt*rinp3;

%% cell morphology
nm(1) = 1; % number of morph #1 (cell body) per cell
ng(1) = 1; % number of compartments per morph

% NOTE: ng for dendrites constrained by k (next section), must be >= 2*k+1
nm(2) = 2; % number of morph #2 (straight dendrites) per cell $$$$$
ng(2) = 20; % number of compartments per morph

dx = 1; % center-to-center spacing of grid compartments (um)

dB = 1/12; % scaling factor for Ca fluxes into cell body

%% Constants for dendritic diffusion computations
k = 2; % Order of Laplacian computation

% for straight dendrites:
L2 = lap(k, ng(2), dx); % matrix to compute Laplacian (straight dendrites)
LN = robinBC(k, ng(2), dx, 0, 1); % get kth-order Neumann BC constants
BN2 = -LN(1,2:end)/LN(1,1); % matrix to find boundary value for 0 slope
LN2 = LN(end,:);  % matrix to compute slope @ proximal end

clear LN;

%% Define cell arrays to hold state & input data & positions

STATES = cell(ncells,nmorphs,nstates);
% index 1: cell number;
% index 2: morphological class:
%   1 => cell body; 2 => straight dendrites;
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
            if m == 1 % for cell body (morph 1)
                STATES{n,m,l} = IC(l)*ones(ng(m),nm(m));
            else % for dendrites (morphs 2 & 3)
                if (l==1 || l==2) % diffusive states have ng+2 grid points
                    STATES{n,m,l} = IC(l)*ones(ng(m)+2,nm(m));
                else % other states have ng grid points
                    STATES{n,m,l} = IC(l)*ones(ng(m),nm(m));
                end
            end
        end
    end
end

% Storing Ca throughout the cells
STATES_store  = cell(itend,ncells,nmorphs);
for ii = 1:itend
    for n = 1:ncells
        for m = 1:nmorphs
            if m == 1 % for cell body (morph 1)
                STATES_store{ii,n,m} = ones(1,1);
            else
                STATES_store{ii,n,m} = ones(ng(m)+2,nm(m));
            end
        end
    end
end

INPUTS = cell(ncells,nmorphs);
% index 1: cell number;
% index 2: morphological class:
%   1 => cell body; 2 => straight dendrites; (3 => branched dendrites);
for n = 1:ncells
    for m = 1:nmorphs
        if m == 1 % for cell body (morph 1)
            INPUTS{n,m} = 0;
        else % for dendrites (morph 2)
            INPUTS{n,m} = zeros(ng(m),nm(m));
        end
    end
end

NGJ = 11; % compartment number for gap junction (both cells)
NGJG = NGJ+1; % grid number for gap junction (used in STATE array)

INTERCONX = [1,2,2,NGJ,2,2,1,NGJ]; % interconnection data, to contain:
% cell #, morph #, segment (column) #, compartment #, for both cells

% POSITIONS indices: cell#, morph#, x or y
POSITIONS = cell(ncells,nmorphs,2);
LD = ng(2)*dx;
FS = [ 0:dx:LD ]';
DL = 2*(LD-(NGJ-1)*dx);
% dendrites, cell 1
POSITIONS{1,2,1} = [ FS, 2*LD-FS ];
POSITIONS{1,2,2} = zeros(ng(2)+1,2);
% body, cell 1
POSITIONS{1,1,1} = LD;
POSITIONS{1,1,2} = 0;
% dendrites, cell 2
POSITIONS{2,2,1} = POSITIONS{1,2,1}+DL;
POSITIONS{2,2,2} = zeros(ng(2)+1,2);
% body, cell 2
POSITIONS{2,1,1} = LD+DL;
POSITIONS{2,1,2} = 0;

% grids for plots
xgridLl = [ POSITIONS{1,2,1}(1:end-1,1) ]';
xgridLr = [ POSITIONS{1,2,1}(1:end-1,2) ]';
xgridRl = [ POSITIONS{2,2,1}(1:end-1,1) ]';
xgridRr = [ POSITIONS{2,2,1}(1:end-1,2) ]';

%% Loop on time
count=1;
for t = 0:dt:tend % loop on time
    
    % set up figure
    figure(1)
    
    plot(xgridLl,STATES{1,2,1}(2:end-1,1),'b')
    hold on
    %plot(xgridLl,STATES{1,2,2}(2:end-1,1),'r')
    plot(POSITIONS{1,1,1},STATES{1,1,1},'b^')
%     plot(POSITIONS{1,1,1},STATES{1,1,2}, 'r^')
    plot(xgridLr,STATES{1,2,1}(2:end-1,2),'b')
    %plot(xgridLr,STATES{1,2,2}(2:end-1,2),'r')
       
    plot(xgridRl,STATES{2,2,1}(2:end-1,1),'Color',[0.4,0.4,1])
    %plot(xgridRl,STATES{2,2,2}(2:end-1,1),'Color',[1,0.4,0.4])
    plot(POSITIONS{2,1,1},STATES{2,1,1},...
        'Marker','^','MarkerEdgeColor',[0.4,0.4,1])
%     plot(POSITIONS{2,1,1},STATES{2,1,2},...
%         'Marker','^','MarkerEdgeColor',[1,0.4,0.4])
    plot(xgridRr,STATES{2,2,1}(2:end-1,2),'Color',[0.4,0.4,1])
    %plot(xgridRr,STATES{2,2,2}(2:end-1,2),'Color',[1,0.4,0.4])

    
    grid on
    xlabel('Distance (uM)');
    ylabel('Concentration (uM)');
    set(gca,'fontsize',15);
    xlim([0 POSITIONS{2,2,1}(1,2)]);
    ylim([0 5]);
    hold off

    % ------------------------------------------------
    % Interconnection between cells [Gap junctions]
    
    % transmission of InP3 through Gap junctions
    IC1 = STATES{1,2,2}(NGJG,2); % [Ca] in compartment in cell 1
    IC2 = STATES{2,2,2}(NGJG,1); % [Ca] in compartment in cell 2
    flow = rinp3 * (IC1 - IC2);
    STATES{1,2,2}(NGJG,2) = STATES{1,2,2}(NGJG,2) - flow;
    STATES{2,2,2}(NGJG,1) = STATES{2,2,2}(NGJG,1) + flow;

    % ------------------------------------------------
    
    % set up inputs: G-protein pulse in 1st 4 compartments of dendrite 1
    if t<0.1 % 100ms long
        INPUTS{1,2}(1,1) = kGI;
%         INPUTS{1,2}(2,1) = kGI;
%         INPUTS{1,2}(3,1) = kGI;
%         INPUTS{1,2}(4,1) = kGI;
    else
        INPUTS{1,2}(1,1) = 0;
%         INPUTS{1,2}(2,1) = 0;
%         INPUTS{1,2}(3,1) = 0;
%         INPUTS{1,2}(4,1) = 0;
    end

    % Following: loops for state updates in individual cells & morphs
    for n = 1:ncells % loop on cells
        for m = 1:nmorphs % loop on morphs
            
            %%%% Implement diffusion of Ca and InP3
            
            % retrieve current calcium and InP3 concentrations
            Ca = STATES{n,m,1};
            InP3 = STATES{n,m,2};

            if m == 1 % diffusion flow for cell body (morph 1)
                % diffusion:
                % grad at proximal ends, straight dendrites
                GCa = LN2*STATES{n,2,1};
                GInP3 = LN2*STATES{n,2,2};
                % Sum fluxes (proportional to grads)
                Ca = Ca - DCa*dB*sum(GCa');
                InP3 = InP3 - DInP3*dB*sum(GInP3');
                % For use in computing boundary conditions
                CaB = Ca; % for use in computing boundary conditions
                InP3B = InP3;  % for use in computing boundary conditions
                % for use in computing internal states
                CaM = Ca;
                InP3M = InP3;
                               
            elseif m == 2 % for straight dendrites
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
            if m==1
                STATES{n,m,1} = CaM;
                STATES{n,m,2} = InP3M;
            else
                STATES{n,m,1} = [Ca(1,:);CaM;Ca(end,:)];
                STATES{n,m,2} = [InP3(1,:);InP3M;InP3(end,:)];
            end
            STATES{n,m,3} = R;
            STATES{n,m,4} = O;
            STATES{n,m,5} = A;
            STATES{n,m,6} = I2;
            STATES{n,m,7} = AA;
            STATES{n,m,8} = CaCalB;
            STATES{n,m,9} = PKCv;
            
            
            STATES_store{count,n,m} = STATES{n,m,1};
            
        end
    end
    pause(0.0002)
    count = count+1;
end

save('Ca_store.mat','STATES_store','-v7.3');
save('Data_files.mat','POSITIONS','INTERCONX');
