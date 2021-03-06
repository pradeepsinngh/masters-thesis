% Plot Ca

dt = 5E-4;   % time step for temporal integration
nmorphs = 3;

% load data files
load('Ca_store.mat');
load('Data_files');
itend = length(STATES_store);
ncells = XE*YE;

% setup scale
XS = [];
YS = [];

for n=1:ncells
    for m=2:nmorphs
        XS = [XS;max(max(POSITIONS{n,m,1})')];
        YS = [YS;max(max(POSITIONS{n,m,2})')];
    end
end

XMAX = ceil(max(XS));
YMAX = ceil(max(YS));
XS = [];
YS = [];

for n=1:ncells
    for m=2:nmorphs
        XS = [XS;max(max(POSITIONS{n,m,1})')];
        YS = [YS;max(max(POSITIONS{n,m,2})')];
    end
end

XMIN = floor(min(XS));
YMIN = floor(min(YS));


figure(3);
%axis square
%axis([XMIN XMAX YMIN YMAX]);
axis([0 XMAX -5 YMAX]);
hold on

tic % timer to cal elapsed time

for ii=1:itend
    for n=1:ncells
        for m=1:nmorphs
            
            % following plots astrocyte skeleton, colored according to [Ca]
            % COLOR contains a 64 x 3 colormap.
            Camax = 9;
            COLOR = colormap;
            
            Ca = STATES_store{ii,n,m};
            Ca = transpose(Ca);
            C = round(Ca*63/Camax)+1; % indices into colormap   
            X = POSITIONS{n,m,1}; % x-coordinates of compartment boundaries
            Y = POSITIONS{n,m,2}; % y-coordinates of compartment boundaries
            
            if ~mod(ii,200)
                if m==1 % for cell body
                                        
                    CB = C;                               
                    XB = X;
                    YB = Y;  
            
                    plot(XB,YB, 'o', 'MarkerSize',1, ...
                        'MarkerEdgeColor',COLOR(CB,:), ...
                        'MarkerFaceColor',COLOR(CB,:));
                    hold on;
                else % for dendrites                  
                    
                    [nr,nc] = size(X);
                    for kk = 1:nc-1
                        for jj = 1:nr
                            plot([X(jj,kk), X(jj,kk+1)], ...
                                [Y(jj,kk), Y(jj,kk+1)], ...
                                'Color', COLOR(C(jj,kk+1),:), ...
                                'LineWidth', 3);
                        end
                    end
                end
                pause(.00002)
            end
        end
    end
end

toc % end timer here