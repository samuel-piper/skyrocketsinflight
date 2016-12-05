function [Avarmax, Rvarmax, Pvarmax, Cvarmax] = Project_2_Sensitivity()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Purpose: To vary the performance dependent parameters in order to
% determine what combinations will result in an accurate landing
%
% Inputs:   None
%
% Outputs:  Avarmax - The maximum values attained while varying launch
%               angle from 0.75 to 1.5 of original angle.
%           Rvarmax - The maximum values attained while varying the volume
%               ratio from 0.75 to 1.5 of original value.
%           Pvarmax - The maximum values attained while varying launch
%               pressure from 0.75 to 1.5 of original value.
%           Cvarmax - The maximum values attained while varying the
%               coefficient of drag from 0.75 to 1.5 of original value.
%         FOR ALL - First column is max height, second column is max
%                       distance
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

% Initialize test matrices
Avarmax = zeros(4,2);
Cvarmax = zeros(4,2);
Pvarmax = zeros(4,2);
Rvarmax = zeros(4,2);


% Reset Variables
b = 1;
c = 1;
d = 1;
i = 1;
figure(2)
hold on
% Vary the performance dependent variables

for a = 0.75:0.25:1.5
    %***** MUST BE UPDATED IN FUNCTION vdp1 AS WELL*****
    Langle = a * pi / 4; % Launch angle (rad)
    pAi = 600000 * b; % Initial pressure of the air (Pa)
    vRatio = 0.8 * c; % Initial volume ratio of the water and air
    Cd = 0.5 * d; % Drag coefficient for the rocket(changes with each model)
    %****************************************************
    R = 287; % Ideal gas constant (J/kg/K)
    
    % Determine initial conditions
    
    TAi = 300; % Initial temperature of the air in the bottle (K)
    %vAi = vB - mWi / rhow; % Initial volume of the air (m^3)
    mB = 0.07; % Mass of the empty bottle, vanes and nose cone.
    mWi = 2 * vRatio; % Initial mass of the water (kg)
    z0 = 0.1; % Initial vertical position (m)
    x0 = 0; % Initial horizontal position (m)
    V0 = 0; % Initial launch velocity (m/s)
    vAi = 0.002 * vRatio; % Initial volume of the air (m^3)
    mAi = pAi * vAi / (R * TAi); % Initial mass of the air (needed for mRi)
    mRi = mB + mAi + mWi; % Initial mass of the entire assembled rocket
    
    % Run differential equations
    
    [~, H] = ode45(@vdp1,[0 25],[z0 x0 mRi V0 Langle vAi mWi mAi]);
    
    % Determine if the results are within tolerance of the target
    
    if max(H(:,2)) > 84 && max(H(:,2)) < 86
        wins(winner) = max(H(:,2));
        winconditions(winner,:) = [Langle, pAi, vRatio, Cd];
        winner = winner + 1;
    end
    
    %Plot the results
    
    plot(H(:,2),H(:,1))
    title('Launch Angle Varied Flight Path')
    xlabel('Horizontal Distance(m)')
    ylabel('Height(m)')
    grid on
    
    %Record results
    Avarmax(i,:) = max([H(:,1),H(:,2)]);
    i = i + 1;
end


% Reset Variables
hold off
a = 1;
c = 1;
d = 1;
i = 1;
figure(3)
hold on

for b = 0.75:0.25:1.5
    %***** MUST BE UPDATED IN FUNCTION vdp1 AS WELL*****
    Langle = a * pi / 4; % Launch angle (rad)
    pAi = 600000 * b; % Initial pressure of the air (Pa)
    vRatio = 0.8 * c; % Initial volume ratio of the water and air
    Cd = 0.5 * d; % Drag coefficient for the rocket(changes with each model)
    %****************************************************
    R = 287; % Ideal gas constant (J/kg/K)
    
    % Determine initial conditions
    
    TAi = 300; % Initial temperature of the air in the bottle (K)
    %vAi = vB - mWi / rhow; % Initial volume of the air (m^3)
    mB = 0.07; % Mass of the empty bottle, vanes and nose cone.
    mWi = 2 * vRatio; % Initial mass of the water (kg)
    z0 = 0.1; % Initial vertical position (m)
    x0 = 0; % Initial horizontal position (m)
    V0 = 0; % Initial launch velocity (m/s)
    vAi = 0.002 * vRatio; % Initial volume of the air (m^3)
    mAi = pAi * vAi / (R * TAi); % Initial mass of the air (needed for mRi)
    mRi = mB + mAi + mWi; % Initial mass of the entire assembled rocket
    
    % Run differential equations
    
    [~, H] = ode45(@vdp1,[0 25],[z0 x0 mRi V0 Langle vAi mWi mAi]);
    
    % Determine if the results are within tolerance of the target
    
    if max(H(:,2)) > 84 && max(H(:,2)) < 86
        wins(winner) = max(H(:,2));
        winconditions(winner,:) = [Langle, pAi, vRatio, Cd];
        winner = winner + 1;
    end
    
    %Plot the results
    
    plot(H(:,2),H(:,1))
    title('Initial Air Pressure Varied Flight Path')
    xlabel('Horizontal Distance(m)')
    ylabel('Height(m)')
    grid on
    
    %Record results
    Pvarmax(i,:) = max([H(:,1),H(:,2)]);
    i = i + 1;
end


% Reset Variables
hold off
a = 1;
b = 1;
d = 1;
i = 1;
figure(4)
hold on

for c = 0.75:0.25:1.5
    %***** MUST BE UPDATED IN FUNCTION vdp1 AS WELL*****
    Langle = a * pi / 4; % Launch angle (rad)
    pAi = 600000 * b; % Initial pressure of the air (Pa)
    vRatio = 0.8 * c; % Initial volume ratio of the water and air
    Cd = 0.5 * d; % Drag coefficient for the rocket(changes with each model)
    %****************************************************
    R = 287; % Ideal gas constant (J/kg/K)
    
    % Determine initial conditions
    
    TAi = 300; % Initial temperature of the air in the bottle (K)
    %vAi = vB - mWi / rhow; % Initial volume of the air (m^3)
    mB = 0.07; % Mass of the empty bottle, vanes and nose cone.
    mWi = 2 * vRatio; % Initial mass of the water (kg)
    z0 = 0.1; % Initial vertical position (m)
    x0 = 0; % Initial horizontal position (m)
    V0 = 0; % Initial launch velocity (m/s)
    vAi = 0.002 * vRatio; % Initial volume of the air (m^3)
    mAi = pAi * vAi / (R * TAi); % Initial mass of the air (needed for mRi)
    mRi = mB + mAi + mWi; % Initial mass of the entire assembled rocket
    
    % Run differential equations
    
    [~, H] = ode45(@vdp1,[0 25],[z0 x0 mRi V0 Langle vAi mWi mAi]);
    
    % Determine if the results are within tolerance of the target
    
    if max(H(:,2)) > 84 && max(H(:,2)) < 86
        wins(winner) = max(H(:,2));
        winconditions(winner,:) = [Langle, pAi, vRatio, Cd];
        winner = winner + 1;
    end
    
    %Plot the results
    
    plot(H(:,2),H(:,1))
    title('Volume Ratio Varied Flight Path')
    xlabel('Horizontal Distance(m)')
    ylabel('Height(m)')
    grid on
    
    %Record results
    Rvarmax(i,:) = max([H(:,1),H(:,2)]);
    i = i + 1;
end


% Reset Variables
hold off
a = 1;
b = 1;
c = 1;
i = 1;
figure(5)
hold on

for d = 0.75:0.25:1.5
    %***** MUST BE UPDATED IN FUNCTION vdp1 AS WELL*****
    Langle = a * pi / 4; % Launch angle (rad)
    pAi = 400000 * b; % Initial pressure of the air (Pa)
    vRatio = 0.8 * c; % Initial volume ratio of the water and air
    Cd = 0.5 * d; % Drag coefficient for the rocket(changes with each model)
    %****************************************************
    R = 287; % Ideal gas constant (J/kg/K)
    
    % Determine initial conditions
    
    TAi = 300; % Initial temperature of the air in the bottle (K)
    %vAi = vB - mWi / rhow; % Initial volume of the air (m^3)
    mB = 0.07; % Mass of the empty bottle, vanes and nose cone.
    mWi = 2 * vRatio; % Initial mass of the water (kg)
    z0 = 0.1; % Initial vertical position (m)
    x0 = 0; % Initial horizontal position (m)
    V0 = 0; % Initial launch velocity (m/s)
    vAi = 0.002 * vRatio; % Initial volume of the air (m^3)
    mAi = pAi * vAi / (R * TAi); % Initial mass of the air (needed for mRi)
    mRi = mB + mAi + mWi; % Initial mass of the entire assembled rocket
    
    % Run differential equations
    [~, H] = ode45(@vdp1,[0 25],[z0 x0 mRi V0 Langle vAi mWi mAi]);
    
    % Determine if the results are within tolerance of the target
    
    if max(H(:,2)) > 84 && max(H(:,2)) < 86
        wins(winner) = max(H(:,2));
        winconditions(winner,:) = [Langle, pAi, vRatio, Cd];
        winner = winner + 1;
    end
    
    %Plot the results
    
    plot(H(:,2),H(:,1))
    title('Drag Coefficient Varied Flight Path')
    xlabel('Horizontal Distance(m)')
    ylabel('Height(m)')
    grid on
    
    %Record results
    Cvarmax(i,:) = max([H(:,1),H(:,2)]);
    i = i + 1;
end
legend('First Launch (0.75 value)','Second Launch (1 value)','Third Launch (1.25 value)','Fourth Launch (1.5 value)')
hold off
end