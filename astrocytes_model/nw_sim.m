% Compute and store states for a cell with multiple dendrites 

addpath('/home/singh/Matlab/mole-master/mole_MATLAB') % path of mole library

dt = 5E-4;   % time step for temporal integration
tend = 5;    % end time ( t(start) = 0)

itend = floor(tend/dt)+1;

%load('IC.mat'); % load initial conditions
% load rate constants / functions for InP3 receptors
load('SD_rates.mat');
load('Data_files.mat');
NN = size(POSITIONS);
ncells = NN(1); 


nmorphs = 3; % number of different component morphologies in cell $$$$$
nstates = 9; % number of relevant state variables

%% ----------------------------------------------------
% constants related to stimulus:

rmax = 2;           % radius of influence 
rmax2 = rmax^2;
tstim = 1/8 * tend; % terminate stimulus at halway through simulation
vel = YE/tend;      % stimulus speed
Xstim = 2;          % x-location of stimulus in 3x3 array 

%% ----------------------------------------------------
% Rate & related constants:

% rate of InP3 production induced by glutamatergic input
% kGI = 100;

DCa = 240;   % Ca diffusion coefficient (units um^2.s^-1)
fCa = 0.3;   % fraction to reduce DCa due to intracellular crowding
DInP3 = 300; % InP3 diffusion coefficient (units um^2.s^-1)
%%%% fInP3 set to zero => no InP3 diffusion %%%%
fInP3 = 0.7; % fraction to reduce DInP3 due to intracellular crowding

% for [Ca]-dependent production of InP3
% kd1f = 6;
% for conversion of InP3 into AA (and thence back into PIP2)
k1b = 2.5;
% for inhibition of InP3 production by PKC
% ki2 = 0.0943; % units uM^-1 (multiplies [PKC*])  (k2 in Kang-Othmer)
% for activation of PKC
k4f = 0.6;    % units s^-1.uM^-1 (multiplies [Ca])    (k4 in Kang-Othmer)
% for deactivation of PKC*    (k5 in Kang-Othmer)
k4b = 0.5;    % units s^-1 
% Calcium infux rate constant for open InP3R channels
% kCa = 500;    % units uM.s^-1 (multiplies effective open-channel fraction)
% Calcium influx rate for open ARC channels
% kaa = 5;      % units uM.s^-1
% Calcium 1st-order pump constants
kp1 = 75;     % units uM.s^-1
kpc1 = 1.8;   % units uM
% Calcium 2nd-order pump constants
kp2 = 12;     % units  units uM.s^-1
kpc2 = 0.10;  % units uM
kpc2sq = kpc2^2;

% rate constant at interconnections
% rinp3 = 500;

% leakage rate
Klk = 1;      % units uM.s^-1

% rates for calcium buffering (with CalB0=40, gives 5/7 buffered @ equil.)
% kcbf = 0.7;   % units s^-1.uM^-1 (multiplies [Ca] & [CalB])
% kcbb = 10;    % units s^-1 (multiplies [CaCalB])

% relative dependence of open probablity on O and A states
fO = 0.1;
fA = 0.9;

PKC0 = 1.0;   % available PKC concentration (uM)   (KT in Kang-Othmer)
CalB0 = 40;   % available Calbindin conc. (uM)   (btot in Breit-Queisser)



%% ------------------------------------------------------------

% local effect
kGI = 400;
kaa = 6;
rinp3 = 500; %5000

% global effect
kd1f = 6; %6
kCa = 600; % 600

kcbf = 0.7;
kcbb = 10;
ki2 = 0.0943;

%% ------------------------------------------------------------
% Scale all rate constants by dt for purposes of temporal integration

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
rinp3 = dt*rinp3;

%% -----------------------------------------------------------------

% cell morphology
nm(1) = 1; % number of morph #1 (cell body) per cell
ng(1) = 1; % number of compartments per morph

% NOTE: ng for dendrites constrained by k (next section), must be >= 2*k+1
NN = size(POSITIONS{1,2,1});
NSS = NN(1);
nm(2) = NN(1); %8 number of morph #2 (straight dendrites) per cell $$$$$
ng(2) = NN(2)-1; % number of compartments per morph

NN = size(POSITIONS{1,3,1});
NBS = NN(1); % number of branched dendrites per cell
nm(3) = NN(1); % 4 number of segments per branched dendrite
ng(3) = NN(2)-1; % number of morph #3 (branched dendrite segments) per cell
ns = 3; % number of dendritic compartments per branched dendrite

dx = 1; % center-to-center spacing of grid compartments (um)

dB = 1/12; % scaling factor for Ca fluxes into cell body

%% ----------------------------------------------------

% Constants for dendritic diffusion computations
k = 2; % Order of Laplacian computation

% for straight dendrites:
L2 = lap(k, ng(2), dx); % matrix to compute Laplacian (straight dendrites)
LN = robinBC(k, ng(2), dx, 0, 1); % get kth-order Neumann BC constants
BN2 = -LN(1,2:end)/LN(1,1); % matrix to find boundary value for 0 slope
LN2 = LN(end,:);  % matrix to compute slope @ proximal end

% for branched dendritic segments:
L3 = lap(k, ng(3), dx); % matrix to compute Laplacian (branched segments)
LN = robinBC(k, ng(3), dx, 0, 1); % get kth-order Neumann BC constants
BN3 = -LN(1,2:end)/LN(1,1); % matrix to find boundary value for 0 slope
LN3 = LN(end,:);  % matrix to compute slope @ proximal end

clear LN;

%% ----------------------------------------------------
% Define cell array to hold state data, and fill with initial conditions

STATES = cell(ncells,nmorphs,nstates);

% itend 1: cell number;
% itend 2: morphological class:
%   1 => cell body; 2 => straight dendrites; 3 => branched dendrites;
% itend 3: states:
%   1 => Ca; 2 => InP3;
%   3 => R; 4 => A; 5 => O; 6 => I2;
%   7 => AA; 8 => CaCalB; 9 => PKCv;

% initial conditions of states for all compartments
IC = [0, 0, 1, 0, 0, 0, 0, 0, 0];

for n = 1:ncells
    for m = 1:nmorphs
        for l = 1:nstates
            if m == 1 % for cell body (morph 1)
                STATES{n,m,l} = IC(l)*ones(ng(m),nm(m));
            else % for dendrites (morphs 2 & 3)
                if (l == 1 || l == 2) % diffusive states (Ca) & InP3 have ng+2 grid points
                    STATES{n,m,l} = IC(l)*ones(ng(m)+2,nm(m));
                else % other states have ng grid points
                    STATES{n,m,l} = IC(l)*ones(ng(m),nm(m));
                end
            end
        end
    end
end

%% ----------------------------------------------------

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

% to hold input (G-protien) inputs
INPUTS = cell(ncells, nmorphs);
% itend 1: cell number
% itend 2: morph class
% 1 => cell body;
% 2 => straight den
% 3 => branched den
for n = 1:ncells
    for m = 1:nmorphs
        if m ==1 % for cell body (morph 1)
            INPUTS{n,m} = 0;
        else
            INPUTS{n,m}= zeros(ng(m), nm(m));
        end
    end
end

% location to recieve (singgle implusive) inputs
% CELL BODY OF MIDDLE CELL
NC = ceil(ncells/2);
NM = 1; % morph number
NG = 1; % grid number
NI = 1; % instance of morph

%% ----------------------------------------------------

% % To store Ca and InP3 history for two compartments.
% 
% % Compartment 1 = Cell 3, SD 5, Comp 1
% % Compartment 2 = Cell 1, BD 2, Comp 3
% 
% % column 1: dt (time step)
% % column 2: Ca (for comp 1)
% % column 3: Ca (for comp 2)
% % column 4: InP3 (for comp 1)
% % column 5: InP3 (for comp 2)
% 
% % store all states for two compartments
% All_STATES = zeros(itend,4,nstates);
% 
% % store AA/ InP3 states for cell 1 compartment 
% AA_STATES = zeros(itend,5);
% 
% InP3_STATES = zeros(itend,2,7);

%% ----------------------------------------------------
% INPUTS = cell(ncells,nmorphs);

% for n = 1:ncells
%     for m = 1:nmorphs
%         if m == 1 % for cell body (morph 1)
%             INPUTS{n,m} = 0;
%         else % for dendrites (morphs 2 & 3)
%             INPUTS{n,m} = zeros(ng(m),nm(m));
%         end
%     end
% end

%% --------------------------------------------------------
% Loop on time

% nsize = size(POSITIONS{1,2,1});
% NSS = nsize(1); % this gives the number of straight dendrites
% nsize = size(POSITIONS{1,3,1});
% NBS = nsize(1); % this gives the number of branched dendritic segments

tic % timer to cal elapsed time


count  = 1;
for t = 0:dt:tend % loop on time
    
    % set  up inputs: G-protien pulse in center cel body'
    if t<0.1 % 100ms long
        INPUTS{NC,NM}(NG,NI) = kGI;
    else
        INPUTS{NC,NM}(NG,NI) = 0;
    end
    
%     Ystim = vel*t; % position of input stimulus
    
    %----------------------------------------------------
    % Interconnection between cells uisng IC.mat
    % load('IC.mat'); % This is done above.
    %
        if isempty(INTERCONX) == 0
            for ii = 1:length(INTERCONX)
    
                % 1st part of interconnection: Cell 1
                cell_1 = INTERCONX(ii,1);
                morph_1 = INTERCONX(ii,2);
                den_1 = INTERCONX(ii,3);
                comp_1 = INTERCONX(ii,4);
    
                % 2nd part of interconnection: Cell 2
                cell_2 = INTERCONX(ii,5);
                morph_2 = INTERCONX(ii,6);
                den_2 = INTERCONX(ii,7);
                comp_2 = INTERCONX(ii,8);
    
                flow = rinp3 * (STATES{cell_1,morph_1,2}(comp_1,den_1) - STATES{cell_2,morph_2,2}(comp_2,den_2));
                STATES{cell_1,morph_1,2}(comp_1,den_1) = STATES{cell_1,morph_1,2}(comp_1,den_1) - flow;
                STATES{cell_2,morph_2,2}(comp_2, den_2) = STATES{cell_2,morph_2,2}(comp_2, den_2) + flow;
            end
        end


   %-----------------------------------------------------
   
%     % Calcium influx from inputs
%     if t < tstim % execute only if stimulus is present
%         for n = 1:ncells % loop on cells

%             Xcell = POSITIONS{n,1,1}; % x of cell center
%             Ycell = POSITIONS{n,1,2}; % y of cell center
%             
%             % following only looks at cells within radius rmax of stimulus
%             if (Xstim-Xcell)^2+(Ystim-Ycell)^2 < rmax2 % stimulus closeness
%                 
%                 for nd = 1:NSS % loop over straight dendrites
%                     Xtip = POSITIONS{n,2,1}(1,1); % x of dendrite tip
%                     Ytip = POSITIONS{n,2,2}(1,1); % y of dendrite tip
%                     
%                     %if stimulus is near this compartment
%                     if (Xstim-Xtip)^2+(Ystim-Ytip)^2 < rmax2
%                         INPUTS{n,2}(1,1) = kGI;
%                     else
%                         INPUTS{n,2}(1,1) = 0;
%                     end
%                 end
%                 
%                 
%                 for nd = 1:NBS % loop over branched dendritic segments
%                     if mod(nd,3) % only execute for terminal segments
%                         Xtip = POSITIONS{n,3,1}(1,1); % x of dendrite tip
%                         Ytip = POSITIONS{n,3,2}(1,1); % y of dendrite tip
%                         
%                         %if stimulus is near this compartment
%                         if (Xstim-Xtip)^2+(Ystim-Ytip)^2 < rmax2
%                             INPUTS{n,3}(1,1) = kGI;
%                         else
%                             INPUTS{n,3}(1,1) = 0;
%                         end
%                     end
%                 end
%             end % end of if clause for stimulus position relative to cell
%         end % end of loop on cells
%     end % end of if clause on t (for stimulus presence)

    

    %----------------------------------------------------
   
    % Following: structure for state updates due to internal processes
    for n = 1:ncells % loop on cells
        for m = 1:nmorphs % loop on morphs
            
            %%%% Implement diffusion of Ca and InP3
            % retrieve current calcium and InP3 concentrations
            Ca = STATES{n,m,1};
            InP3 = STATES{n,m,2};
            
            if m == 1
                CaM = Ca;
                InP3M = InP3;
            else
                CaM = Ca(2:end-1,:);
                InP3M = InP3(2:end-1,:);
            end            

            % retrieve states
            R = STATES{n,m,3};
            O = STATES{n,m,4};
            A = STATES{n,m,5};
            I2 = STATES{n,m,6};
            AA = STATES{n,m,7};
            CaCalB = STATES{n,m,8};
            PKCv = STATES{n,m,9};
            
            %----------------------------------------------------
            %%%% Implement reactions & channel flows
         
            % retrieve inputs
            Gin = INPUTS{n,m};
            
            % InP3 update
            InP3M = InP3M + (Gin + kd1f*CaM)./(1+ki2*PKCv) - k1b*InP3M;
            
            % AA update
            AA = AA + k1b*InP3M - k1b*AA;
            
            % PKCv update
            PKCv = PKCv + k4f*CaM.*(PKC0-PKCv) - k4b*PKCv;
            
            % indices into tables of Ca-dependent rates
            iCa = round(5E2*CaM)+1;

            % InP3 channel states
            DOR = ph2bA(iCa).*O;
            DRO = ph2fA(iCa).*InP3M.*R;
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
            
            CaM(CaM < 0) = 0;
            InP3M(InP3M < 0) = 0;
            
            
             % replace arrays in STATES
            if m==1
                Ca=  CaM;
                InP3 = InP3M;
            else
                STATES{n,m,1} = [Ca(1,:);CaM;Ca(end,:)];
                Ca = STATES{n,m,1};
                STATES{n,m,2} = [InP3(1,:);InP3M;InP3(end,:)];
                InP3 = STATES{n,m,2};
            end


            % Diffusion Process
            % -----------------------------------
            if m == 1 % efflux & diffusion process for cell body (morph 1)
                
                % grad(Ca) at proximal ends, straight dendrites
                GCa = LN2*STATES{n,2,1};
                GInP3 = LN2*STATES{n,2,2};
                
                % grad(Ca) at proximal ends, branched dendrites
                GCa = [GCa, LN3*STATES{n,3,1}(:,ns:ns:end)]; 
                GInP3 = [GInP3, LN3*STATES{n,3,2}(:,ns:ns:end)];
                
                % Sum fluxes (proportional to grads) to increment body Ca
                Ca = Ca - DCa*dB*sum(GCa');
                InP3 = InP3 - DInP3*dB*sum(GInP3');
                
                if Ca<0 % in case roundoff error makes [Ca]<0
                    Ca = 0;
                end
                
                CaB = Ca; % for use in computing boundary conditions                
                InP3B = InP3;
                
            else % efflux & diffusion process for dendrites (morphs 2 & 3)
                if m == 2 % for straight dendrites
                    
                    % Ca diffusion eq
                    Ca = Ca + DCa*L2*Ca;            % lateral diffusion
                    Ca(1,:) = BN2*Ca(2:end,:);      % grad (i.e. flux) = 0, distal ends
                    Ca(end,:) = CaB;                % body concentration, proximal ends
                    
                    % inp3 diffusion eq
                    InP3 = InP3 + DInP3*L2*InP3;   % lateral diffsuion
                    InP3(1,:) = BN2*InP3(2:end,:); % boundary condition (distal end)
                    InP3(end,:) = InP3B;           % boundary condition (open proximal end)
                    
                elseif m == 3 % for branched dendrites
                    
                    % lateral diffusion
                    Ca = Ca + DCa*L3*Ca;
                    InP3 = InP3 + DInP3*L3*Ca;
                    
                    % grad (i.e. flux) = 0, distal ends
                    Ca(1,ns-1:ns:end-1) = BN3*Ca(2:end,ns-1:ns:end-1);
                    Ca(1,ns-2:ns:end-2) = BN3*Ca(2:end,ns-2:ns:end-2);
                    
                    InP3(1,ns-1:ns:end-1) = BN3*InP3(2:end,ns-1:ns:end-1);
                    InP3(1,ns-2:ns:end-2) = BN3*InP3(2:end,ns-2:ns:end-2);
                    
                    % concentration at junctions s.t. sum(fluxes) = 0:
                    Ca(1,ns:ns:end) = 0.5*BN3*( Ca(2:end,ns:ns:end) ...
                        + 0.5*(flip( Ca(1:end-1,ns-1:ns:end-1) ...
                        + Ca(1:end-1,ns-2:ns:end-2) )) );
                    Ca(end,ns-1:ns:end-1) = Ca(1,ns:ns:end);
                    Ca(end,ns-2:ns:end-2) = Ca(1,ns:ns:end);
                    
                    
                    InP3(1,ns:ns:end) = 0.5*BN3*( InP3(2:end,ns:ns:end) ...
                        + 0.5*(flip( InP3(1:end-1,ns-1:ns:end-1) ...
                        + InP3(1:end-1,ns-2:ns:end-2) )) );
                    InP3(end,ns-1:ns:end-1) = InP3(1,ns:ns:end);
                    InP3(end,ns-2:ns:end-2) = InP3(1,ns:ns:end);
                    
                    % concentration = body concentration, proximal ends
                    Ca(end,ns:ns:end) = CaB;
                    InP3(end,ns:ns:end) = InP3B;
                    
                end
            end
            
            
            STATES{n,m,1} = Ca;
            STATES{n,m,2} = InP3;
            STATES_store{count,n,m} = Ca;
            STATES{n,m,3} = R;
            STATES{n,m,4} = O;
            STATES{n,m,5} = A;
            STATES{n,m,6} = I2;
            STATES{n,m,7} = AA;
            STATES{n,m,8} = CaCalB;
            STATES{n,m,9} = PKCv;
            
        end       
    end
    count = count + 1;
end

toc

save('Ca_store.mat','STATES_store','NC','NM','NG','NI','-v7.3')   % save variable in the output.mat file


