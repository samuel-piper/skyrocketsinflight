function [T, H] = Project_2_Test()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Purpose:  Provide proof of concept for the code using the conditions
% outlined in the test case.
%
% Inputs:   None           
%
% Outputs:  T - Times for calls to ode45
%           H - Conditions throughout the test
%
% Assumptions:  None
%
% Created: 11/14/16
% Modified: 11/30/16
% Author: a2f341a79180
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set some required initial values

%***** MUST BE UPDATED IN FUNCTION vdp1 AS WELL*****
global TAi
TAi = 300; % Initial temperature of the air in the bottle (K)
global pAi
pAi = 446090.8; % Initial pressure of the air (Pa)
global vRatio
vRatio = .5; % Initial volume ratio of the water and air
global Cd
Cd = 0.5; % Drag coefficient for the rocket(changes with each model)
%****************************************************
R = 287; % Ideal gas constant (J/kg/K)


% Determine initial conditions

mB = 0.07; % Mass of the empty bottle, vanes and nose cone.
mWi = 2 * vRatio; % Initial mass of the water (kg)
z0 = 0.1; % Initial vertical position (m)
x0 = 0; % Initial horizontal position (m)
V0 = 0; % Initial launch velocity (m/s)
Langle = pi / 4; % Launch angle (rad)
vAi = 0.002 * vRatio; % Initial volume of the air (m^3)
mAi = pAi * vAi / (R * TAi); % Initial mass of the air (needed for mRi)
mRi = mB + mAi + mWi; % Initial mass of the entire assembled rocket

% Run differential equations

[T, H] = ode45(@vdp1,[0 5],[z0 x0 mRi V0 Langle vAi mWi mAi]);

%Plot the results

figure(1)
plot(H(:,2),H(:,1),'g')
title('Rocket Flight Path')
xlabel('Horizontal Distance(m)')
ylabel('Height(m)')
grid on



end