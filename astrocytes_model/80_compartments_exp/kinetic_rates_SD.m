%%%% Sneyd-Dufour InP3 receptor model:
%%%% Compute & store rate constants and functions

%% Model and notation:
% Notes:
% f = forward; b = backward (pos. and neg. digits, resp., used in paper);
% forward <=> *away* from R state
%   1st row: S-D paper notation (digits are *subscripts*);
%   2nd row: my original notation
%   k1  k-1 k2  k-2 k3  k-3 k4  k-4 L1  L3  L5  l2  l-2 l4  l-4 l6  l-6
%   k1f k1b k2f k2b k3f k3b k4f k4b L1  L3  L5  l2f l2b l4f l4b l6f l6b
% ph used for [Ca]-dependent rate functions; r for rate constants
% r1b = r5b = k1b+l2b; called k1l2b in prior scripts (versions v0)
% r3b = k3b;  k3b notation used in prior scripts (versions v0)
%
%            ph1f
%       R   ======  I1
%      ||    r1b
%  ph2f||ph2b
% *InP3||    ph4f       ph5f
%       O   ======  A  ======  I2
%      ||    ph4b       r5b
%  ph3f||r3b
%      ||
%       S
%
% Reduced (simplified) model structure:
%
%       R
%      ||
%  ph2f||ph2b
% *InP3||    ph4f       ph5f
%       O   ======  A  ======  I2
%            ph4b       r5b
%    

%% Fundamental channel constants
k1f = 0.640;   % units s^-1.uM^-1
k1b = 0.040;   % units s^-1
k2f = 37.40;   % units s^-1.uM^-1
k2b = 1.400;   % units s^-1
k3f = 0.110;   % units s^-1 IN PAPER: units = s^1.uM^1! Must be wrong
k3b = 29.80;   % units s^-1
k4f = 4.000;   % units s^-1.uM^-1
k4b = 0.540;   % units s^-1
L1 = 0.120;    % units uM
L3 = 0.025;    % units uM
L5 = 57.40;    % units uM
l2f = 1.700;   % units s^-1
l2b = 0.800;   % units s^-1
l4f = 1.700;   % units s^-1.uM^-1
l4b = 2.500;   % units s^-1.uM^-1
l6f = 4707;    % units s^-1
l6b = 11.40;   % units s^-1

%% Rate constants
r1b = k1b + l2b; % derived rate constant; units s^-1
r3b = k3b; % renamed according to convention; units s^-1
r5b = r1b; % renamed according to convention; units s^-1


%% Concentration-dependent rate functions
% values stored in arrays for table look-up
for jj = 1:15001
    
    Ca = 2E-3*(jj-1); % Ca concentrations range 0 to 10uM
    CaA(jj,1) = Ca; % store values of calcium to allow plotting of functions
    
    % all rate functions except ph2f in units s^-1
    % ph2f in units s^-1.uM^-1
    ph1fA(jj,1) = ( (k1f*L1 + l2f)*Ca ) / ( L1 + (1+L1/L3)*Ca );
    ph2fA(jj,1) = ( k2f*L3 + l4f*Ca ) / ( L3 + (1+L3/L1)*Ca );
    ph2bA(jj,1) = ( k2b + l4b*Ca ) / ( 1 + Ca/L5 );
    ph3fA(jj,1) = ( k3f*L5 ) / ( L5 + Ca );
    ph4fA(jj,1) = ( (k4f*L5 + l6f)*Ca ) / ( L5 + Ca );
    ph4bA(jj,1) = ( L1*( k4b + l6b) ) / ( L1 + Ca );
    ph5fA(jj,1) = ( (k1f*L1 + l2f)*Ca ) / ( L1 + Ca );
    
end

%% Save computed data
save SD_rates.mat r1b r3b r5b ph1fA ph2fA ph2bA ph3fA ph4fA ph4bA ph5fA
