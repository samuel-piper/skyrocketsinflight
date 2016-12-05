function [T, H, w, wc] = Project_2_Var()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Purpose: To vary the performance dependent parameters in order to
% determine what combinations will result in an accurate landing
%
% Inputs:   None           
%
% Outputs:  T - Times for calls to ode45
%           H - Conditions throughout the test
%           w - Landing point for each path taken
%           wc - Conditions that lead to a successful landing
%
% Assumptions:  None
%
% Created: 11/14/16
% Modified: 11/30/16
% Author: a2f341a79180
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Declare needed global values
global pAi
global vRatio
global Cd
global TAi
winner = 1;
wins = zeros(1,100000);
winconditions = zeros(100000,4);


% Vary the performance dependent variables

for a = 0.45:0.05:1.3
    for b = 0.7:0.1:1.3
        for c = 0.7:0.1:1.2
            for d = 0.8:0.2:1.2
                %***** MUST BE UPDATED IN FUNCTION vdp1 AS WELL*****
                Langle = a * pi / 4; % Launch angle (rad)
                pAi = 600000 * b; % Initial pressure of the air (Pa)
                vRatio = 0.8 * c; % Initial volume ratio of the water and air
                Cd = 0.5 * d; % Drag coefficient for the rocket(changes with each model)
                %****************************************************
                R = 287; % Ideal gas constant (J/kg/K)
                
                % Determine initial conditions
                
                TAi = 300; % Initial temperature of the air in the bottle (K)
                mB = 0.07; % Mass of the empty bottle, vanes and nose cone.
                mWi = 2 * vRatio; % Initial mass of the water (kg)
                z0 = 0.1; % Initial vertical position (m)
                x0 = 0; % Initial horizontal position (m)
                V0 = 0; % Initial launch velocity (m/s)
                vAi = 0.002 * vRatio; % Initial volume of the air (m^3)
                mAi = pAi * vAi / (R * TAi); % Initial mass of the air (needed for mRi)
                mRi = mB + mAi + mWi; % Initial mass of the entire assembled rocket
                
                % Run differential equations
                
                [T, H] = ode45(@vdp1,[0 25],[z0 x0 mRi V0 Langle vAi mWi mAi]);
                
                % Determine if the results are within tolerance of the target
                
                if max(H(:,2)) > 84.9 && max(H(:,2)) < 85.1
                    wins(winner) = max(H(:,2));
                    winconditions(winner,:) = [Langle, pAi, vRatio, Cd];
                    winner = winner + 1;
                end                
            end
        end
    end
end

% Condition resultant matrices

sizewins = size(wins);
for i = 100:sizewins(2)
    if wins(i) == 0
        w = wins(1:i);
        wc = winconditions(1:i,:);
        break
    end
end
end