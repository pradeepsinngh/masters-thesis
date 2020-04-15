% Compute wave velocity in a wavefront.

dt = 5E-4;   % time step for temporal integration
nmorphs = 3;

% load data files
load('Ca_store.mat');
load('Data_files');
itend = length(STATES_store);
ncells = XE*YE;


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

% Cell array for data to compute speed
% First index covers cells in facilitatory network
% Second index (for each cell):
% = 1: distance between wave origin and current farthest peak compartment
% = 2: time stamp for current peak
% = 3: distance between wave origin and previous farthest peak compartment
% = 4: time stamp for previous peak
Peaks = -ones(ncells, 4);

% Ca threshold: NOTE THIS WILL BE PARAMETER-DEPENDENT
Ca_threshold = 1.5;

% Origin (stimulus location) x and y-cordinates.
X0 = POSITIONS{NC,NM,1}(NI,NG);
Y0 = POSITIONS{NC,NM,2}(NI,NG);

Ca0 = 0; % Ca conc. at origin
CaD = zeros(ncells, nm(2)+2*nm(3)/3); % ca conc. at den tips.
nm0 = nm(2); % final index into CaD for SD tips
nm1 = nm0+1; % starting index into CaD for branched dendrite tips
nm2 = nm0+2*nm(3)/3; % final index into CaD for branched dendrite tips
m3tipF = ~( ~mod(1:nm(3),3)); % tag indices of branched den tips
m3tipB = 1:nm(3);
m3tipB = m3tipB(m3tipF); % indices into POSITIONS of branched dendrite tips
% distances and elapsed times when farthest tip in cell 'lights up'
farPeaks = zeros(ncells,2);

%%  loop on time
for it = 1:itend % time index
    timeStamp = dt * (it-1);
    
    Ca0p = Ca0;
    Ca0 = STATES_store{it, NC, NM}(NG,NI); % current Ca @ origin
    if (Ca0 >= Ca_threshold && Ca0p < Ca_threshold)
        t0 = timeStamp;
    end
    
    CaDp = CaD;
    for n = 1:ncells % loop on cells in network
        
        % load current Ca at den tips
        for m=2:nmorphs % loop on den for each cell
            if m==2 % for SD (morph 2)
                CaD(n,1:nm0) = STATES_store{it,n,m}(1,:);
            else % for BD (morph 3)
                CaD(n,nm1:nm2) = STATES_store{it,n,m}(1,m3tipF);
            end
        end
        
        
        for k=1:nm0 % for SD
            % look t see if Ca crosses threshold in ea tip
            if (CaD(n,k) >= Ca_threshold && CaDp(n,k) < Ca_threshold)
                % if so, compute distance to wave origin
                X = POSITIONS{n,2,1}(k,1)- X0;
                Y = POSITIONS{n,2,2}(k,1) - Y0;
                D = sqrt(X^2 + Y^2);
                if D > farPeaks(n,1)
                    % if this is the most distant peak so far
                    farPeaks(n,1) = D;
                    farPeaks(n,2) = timeStamp - t0;
                end
            end
        end
        
        for k=nm1:nm2 % for SD
            % look t see if Ca crosses threshold in ea tip
            if (CaD(n,k) >= Ca_threshold && CaDp(n,k) < Ca_threshold)
                % if so, compute distance to wave origin
                k0 = m3tipB(k-nm0);
                X = POSITIONS{n,3,1}(k0,1)- X0;
                Y = POSITIONS{n,3,2}(k0,1) - Y0;
                D = sqrt(X^2 + Y^2);
                if D > farPeaks(n,1)
                    % if this is the most distant peak so far
                    farPeaks(n,1) = D;
                    farPeaks(n,2) = timeStamp - t0;
                end
            end
        end
          
    end % end loop over cells

end % end time loop

%% Mean Speed
speeds = farPeaks(:,1)./farPeaks(:,2);
NC = ceil(ncells/2);
speed = mean([speeds(1:NC-1);speeds(NC+1:end)]);