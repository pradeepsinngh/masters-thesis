%%%% GenerateAstrocytesHexRand.m %%%%

% Script generates a hexagonal array of 'astrocytes' 
% with random variations in position and cell orientation,
% and inteconnects them quasi-randomly, as specified.
% Interconnections are specified in terms of which 
% dendritic compartments are connected.
% It also calculates and writes out the spatial coordinates of
% the boundary points of the dendritic compartments

% astrocyte array size parameters
XE = 3; % X extent of array (number of columns of cells)
YE = 3; % Y extent of array (number of rows of cells)
NG = XE*YE; % total number of astrocytes;

% astrocyte morphology parameters: model-dependent!
LD = 10; % physical length of straight dendrite
halfLD = LD/2;
twoLDsqr = (2*LD)^2;
ND = 12; % total # basal 'dendrites' per cell
NB = 4; % # of these which are branched (MUST BE EVEN DIVISOR OF ND)
NF = ND/NB; % 
NS = ND-NB; % # which are straight
threeNB = 3*NB;

% astrocyte spacing parameters
RD = 0.9; % ratio, dendritic length to inter-cell spacing
DY = LD/RD; % mean vertical separation (intercellular distance), hex array
DX = DY * sqrt(6)/2; % mean horizontal separation, hex array
STD = 0.12*DY; % standard deviation of variations in cell positions
%STD = 0.0; % standard deviation of variations in cell positions

NISmax = 2; % max number of interconnections per straight dendrite
NIBmax = NISmax/2; % % max number of interconnections per branched segment
NC = 10; % number of compartments (elements) per length LD of dendrite
halfNC = NC/2;
eps = 1E-6; % bias parameter for assigning compartments in forked dendrites

nmorph = 3; % 3 morphs = cell body, straight dendrites, branched dendrites
% cell body & compartment boundary coordinates
POSITIONS = cell(XE*YE, nmorph, 2); % indices: cell#, morph#, x or y
INTERCONX = []; % interconnection data, to contain:
% cell #, morph #, segment (column) #, compartment #, for both cells

rng(7); % seed random number generator, for reproduceability

% prototype geometric model of astrocyte, centere at origin:
% a cell with ND basal 'dendrites', of length LD, every ND/NB one forked

FS = (LD/NC).*(NC:-1:0)'; % compartment partition, straight dendrite
FB = (halfLD/halfNC).*(halfNC:-1:0)'; % compartment partition, branched den

XCS = []; % x-coordinates, straight dendrite
YCS = []; % y-coordinates, straight dendrite
XCB = []; % x-coordinates, branched dendrite
YCB = []; % y-coordinates, branched dendrite

for jj=1:ND
    th = 2*pi*(jj-1)/ND;
    if mod(jj,NF) % straight dendrites
        XCS = [ XCS, cos(th)*FS ];
        YCS = [ YCS, sin(th)*FS ];
    else  % forked dendrites
        % branch point
        XM = cos(th)*halfLD;
        YM = sin(th)*halfLD;
        % distal branch #1
        XCB = [ XCB, cos(th-pi/12)*FB+XM ];
        YCB = [ YCB, sin(th-pi/12)*FB+YM ];
        % distal branch #2
        XCB = [ XCB, cos(th+pi/12)*FB+XM ];
        YCB = [ YCB, sin(th+pi/12)*FB+YM ];
        % proximal branch
        XCB = [ XCB, cos(th)*FB ];
        YCB = [ YCB, sin(th)*FB ];
    end
end

XC=[];
YC=[];
% generate positions of cell centers (hexagonal grid w/ random variations)
for ii=1:YE % row counter
    Nii = XE*(ii-1);
    for jj=1:XE % column counter
        N = Nii + jj;
        X = jj*DX + random('Normal',0,STD); % X-coordinate
        % Y-coordinates, with vertical offset between alternating columns
        if mod(jj,2)
            Y = ii*DY + random('Normal',0,STD);
        else
            Y = (ii-0.5)*DY + random('Normal',0,STD);
        end
        % put body x and y coordinates into cell array
        XC = [XC;X];
        YC = [YC;Y];
        POSITIONS{N,1,1} = X;
        POSITIONS{N,1,2} = Y;
    end
end

%%%% plot cell centers
figure(5);
scatter(XC,YC, 30,'k','filled');
hold on;

% place dendrites for first astrocyte
th = (pi/2)*rand; % random rotation of dendritic tree
% x and y coordinates of compartment boundaries, straight dendrites
XS = [ cos(th)*XCS-sin(th)*YCS + XC(1) ];
YS = [ sin(th)*XCS+cos(th)*YCS + YC(1) ];
% x and y coordinates of compartment boundaries, branched dendrites
XB = [ cos(th)*XCB-sin(th)*YCB + XC(1) ];
YB = [ sin(th)*XCB+cos(th)*YCB + YC(1) ];
% transpose coordinate arrays before storing in POSITIONS cell array
POSITIONS{1,2,1} = XS';
POSITIONS{1,2,2} = YS';
POSITIONS{1,3,1} = XB';
POSITIONS{1,3,2} = YB';

%%%% plot dendrites
plot(XS, YS, 'k', 'LineWidth', 1.0)
plot(XB, YB, 'k', 'LineWidth', 1.0)

NIS = zeros(NG,NS); % to hold # interconnections per straight dendrite
NIB = zeros(NG,threeNB); % to hold # interconnections per branched dendrite
% specify remaining segments and interconnections
for ii=2:NG
    % place dendrites for subsequent astrocytes
    th = (pi/2)*rand; % random rotation of dendritic tree
    % x and y coordinates of compartment boundaries, straight dendrites
    XS = [ cos(th)*XCS-sin(th)*YCS + XC(ii) ];
    YS = [ sin(th)*XCS+cos(th)*YCS + YC(ii) ];
    % x and y coordinates of compartment boundaries, branched dendrites
    XB = [ cos(th)*XCB-sin(th)*YCB + XC(ii) ];
    YB = [ sin(th)*XCB+cos(th)*YCB + YC(ii) ];
    % transpose coordinate arrays before storing in POSITIONS cell array
    POSITIONS{ii,2,1} = XS';
    POSITIONS{ii,2,2} = YS';
    POSITIONS{ii,3,1} = XB';
    POSITIONS{ii,3,2} = YB';
       
    %%%% plot dendrites
    plot(XS, YS, 'k', 'LineWidth', 1.0)
    plot(XB, YB, 'k', 'LineWidth', 1.0)
    
    % check for overlap with 'dendrites' of previously-defined astrocytes
    for kk=1:ii-1
        % check if cells are close enough for overlap to occur
        if (XC(ii)-XC(kk))^2+(YC(ii)-YC(kk))^2 < twoLDsqr
            
            for jj=1:NS % loop over straight dendrites for cell# ii
                Xjj = POSITIONS{ii,2,1}(jj,:);
                Yjj = POSITIONS{ii,2,2}(jj,:);
                
                for ll=1:NS % loop over straight dendrites for cell# kk
                    Xll = POSITIONS{kk,2,1}(ll,:);
                    Yll = POSITIONS{kk,2,2}(ll,:);
                    % proceed if # interconnx < max for both 'dendrites'
                    if NIS(ii,jj)<NISmax && NIS(kk,ll)<NISmax
                    % determine points @ which dendrites overlap
                        [X, Y] = polyxpoly(Xjj',Yjj', Xll',Yll');
                        if ~isempty(X) % if overlaps are found
                            scatter(X,Y, 8, 'r');
                            
                            % find compartments of ea 'dendrite' @ overlap
                            for nc=1:length(X)
                                % find compartment # on cell #ii dendrite
                                CNjj = ...
                                    ceil( NC*(sqrt( ...
                                    (X(nc)-Xjj(1))^2 + ...
                                    (Y(nc)-Yjj(1))^2) ) / LD );
                                % find compartment # on cell #kk dendrite
                                CNll = ...
                                    ceil( NC*(sqrt( ...
                                    (X(nc)-Xll(1))^2 + ...
                                    (Y(nc)-Yll(1))^2) ) / LD );
                                INTERCONX = [ INTERCONX; ...
                                    ii,2,jj,CNjj, kk,2,ll,CNll];
                            end
                            
                            NIS(ii,jj) = NIS(ii,jj)+length(X);
                            NIS(kk,ll) = NIS(kk,ll)+length(X);

                        end
                    end
                end
                    
                for ll = 1:threeNB % loop over branched dendrites,cell# kk
                    Xll = POSITIONS{kk,3,1}(ll,:);
                    Yll = POSITIONS{kk,3,2}(ll,:);
                    % proceed if # interconnx < max for both 'dendrites'
                    if NIS(ii,jj)<NISmax && NIB(kk,ll)<NIBmax
                        % determine points @ which dendrites overlap
                        [X, Y] = polyxpoly(Xjj',Yjj', Xll',Yll');
                        if ~isempty(X) % if overlaps are found
                            scatter(X,Y, 6, 'r');

                            % find compartments of ea 'dendrite' @ overlap
                            for nc=1:length(X)
                                % find compartment # on cell #ii dendrite
                                CNjj = ...
                                    ceil( halfNC*(sqrt( ...
                                    (X(nc)-Xjj(1))^2 + ...
                                    (Y(nc)-Yjj(1))^2) ) / halfLD );
                                % find compartment # on cell #kk dendrite
                                CNll = ...
                                    ceil( halfNC*(sqrt( ...
                                    (X(nc)-Xll(1))^2 + ...
                                    (Y(nc)-Yll(1))^2) ) / halfLD );
                                INTERCONX = [ INTERCONX; ...
                                    ii,2,jj,CNjj, kk,3,ll,CNll];
                            end
                            
                            NIS(ii,jj) = NIS(ii,jj)+length(X);
                            NIB(kk,ll) = NIB(kk,ll)+length(X);
                            
                        end
                    end
                end
            end
            
            for jj=1:threeNB % loop over branched dendrites for cell #ii
                Xjj = POSITIONS{ii,3,1}(jj,:);
                Yjj = POSITIONS{ii,3,2}(jj,:);
                
                for ll=1:NS % loop over straight dendrites for cell# kk
                    Xll = POSITIONS{kk,2,1}(ll,:);
                    Yll = POSITIONS{kk,2,2}(ll,:);
                    % proceed if # interconnx < max for both 'dendrites'
                    if NIB(ii,jj)<NIBmax && NIS(kk,ll)<NISmax
                    % determine points @ which dendrites overlap
                        [X, Y] = polyxpoly(Xjj',Yjj', Xll',Yll');
                        if ~isempty(X) % if overlaps are found
                            scatter(X,Y, 6, 'r');
                            
                            % find compartments of ea 'dendrite' @ overlap
                            for nc=1:length(X)
                                % find compartment # on cell # ii dendrite
                                CNjj = ...
                                    ceil( NC*(sqrt( ...
                                    (X(nc)-Xjj(1))^2 + ...
                                    (Y(nc)-Yjj(1))^2) ) / LD );
                                % find compartment # on cell # kk dendrite
                                CNll = ...
                                    ceil( NC*(sqrt( ...
                                    (X(nc)-Xll(1))^2 + ...
                                    (Y(nc)-Yll(1))^2) ) / LD );
                                INTERCONX = [ INTERCONX; ...
                                    ii,3,jj,CNjj, kk,2,ll,CNll];
                            end
                            
                            NIB(ii,jj) = NIB(ii,jj)+length(X);
                            NIS(kk,ll) = NIS(kk,ll)+length(X);
                            
                        end
                    end
                end
                    
                for ll = 1:threeNB
                    Xll = POSITIONS{kk,3,1}(ll,:);
                    Yll = POSITIONS{kk,3,2}(ll,:);
                    % proceed if # interconnx < max for both 'dendrites'
                    if NIB(ii,jj)<NIBmax && NIB(kk,ll)<NIBmax
                        % determine points @ which dendrites overlap
                        [X, Y] = polyxpoly(Xjj',Yjj', Xll',Yll');
                        if ~isempty(X) % if overlaps are found
                            scatter(X,Y, 6, 'r');

                            % find compartments of ea 'dendrite' @ overlap
                            for nc=1:length(X)
                                % find compartment # on cell # ii dendrite
                                CNjj = ...
                                    ceil( halfNC*(sqrt( ...
                                    (X(nc)-Xjj(1))^2 + ...
                                    (Y(nc)-Yjj(1))^2) ) / halfLD );
                                % find compartment # on cell # kk dendrite
                                CNll = ...
                                    ceil( halfNC*(sqrt( ...
                                (X(nc)-Xll(1))^2 + ...
                                    (Y(nc)-Yll(1))^2) ) / halfLD );
                                INTERCONX = [ INTERCONX; ...
                                    ii,3,jj,CNjj, kk,3,ll,CNll];
                            end
                            
                            NIB(ii,jj) = NIB(ii,jj)+length(X);
                            NIB(kk,ll) = NIB(kk,ll)+length(X);
                            
                        end
                    end
                end
            end                 
        end    
    end    
end

AXMAX = max(XE*DX,YE*DY)+LD+3*STD;
AXMIN = -0.5*LD-3*STD;
axis square
axis([-3*STD AXMAX -0.5*LD-3*STD AXMAX]);

save('Data_files.mat','XE','YE','POSITIONS','INTERCONX')   % save variable in the output.mat file
