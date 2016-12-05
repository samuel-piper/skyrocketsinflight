function dhdt = vdp1(~, h)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Purpose: Establish differential equations for use as input for call to
% ode45.
%
% Inputs:   t - Bounds from time for call to ode45(not used)
%           h - h(1) is vertical position
%               h(2) is horizontal position
%               h(3) is mass of the entire rocket
%               h(4) is magnitude of the velocity
%               h(5) is tragectory angle
%               h(6) is volume of air in the bottle
%               h(7) is mass of the water
%               h(8) is mass of the air
%               h(9) is a check for landing conditions
%
% Outputs:  dhdt - Differential equations for determining flight path.
%               dhdt(1) is the change in vertical position
%               dhdt(2) is the change in horizontal position
%               dhdt(3) is the change in mass of the entire rocket
%               dhdt(4) is the change in velocity
%               dhdt(5) is the change in tragectory angle
%               dhdt(6) is the change in volume of the air in the bottle
%               dhdt(7) is the change in mass of water
%               dhdt(8) is the change in mass of air
%               dhdt(9) is the trigger for the check on landing conditions
%
% Assumptions:  None
%
% Created: 11/14/16
% Modified: 11/23/16
% Author: a2f341a79180
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize outside parameters

%***** MUST BE UPDATED IN MAIN FUNCTION AS WELL*****
global vRatio % Initial volume ratio of the water and air
mB = 0.07; % Mass of the empty bottle, vanes and nose cone.
global TAi % Initial temperature of the air in the bottle (K)
global pAi % Initial pressure of the air (Pa)
%****************************************************

% Initialize performance dependent variables

% Initial pressure of air (pAi) set above in 'Initialize outside parameters' 
% Volume ratio (vRatio) is set above in 'Initialize outside parameters'
mWi = 2 * vRatio; % Initial mass of the water (kg)
global Cd % Drag coefficient for the rocket(changes with each model)

% Initialize the function variables

dhdt = zeros(8, 1);
z = h(1); % vertical position (m)
x = h(2); % horizontal position (m)
mR = h(3); % mass of the entire rocket (kg)
V = h(4); % velocity (m/s)
theta = h(5); % tragectory angle (rad)
vA = h(6); % volume of air in the bottle (m^3)
mW = h(7); % mass of the water (kg)
mA = h(8); % mass of the air (kg)

% Initialize other variables

g0 = 9.81; % Gravitational acceleration (m/s^2)
R = 287; % Ideal gas constant (J/kg/K)
rhow = 1000; % Density of water (kg/m^3)
vB = 0.002; % Volume of the bottle (m^3)
CircB = 0.32986; % Circumference of the bottle (m)
At = 0.0013854; % Area of the throat of the bottle (m^2)
rhoatm = 0.961; % Density of the surrounding air (kg/m^3)
Ab = pi * ( CircB / (2 * pi ) )^2; % Area of the bottle (m^2)
vAi = 0.002 * vRatio; % Initial volume of the air (m^3)
Patm = 82943.93; % Atmospheric pressure in Boulder, CO (Pa)
Cdd = 0.8; % Discharge coefficient

% Run basic equations

q = rhoatm * V^2 / 2; % Dynamic Pressure
D = q * Cd * Ab; % Drag on the rocket
mAi = pAi * vAi / (R * TAi); % Initial mass of the air 

% Check and correct anomalous conditions
if mW < 0
    mW = 0;
end
if vA > 0.002
    vA = 0.002;
end
if -pi/2 > theta
    theta = -pi/2;
end

% Thermodynamic/mass equations for phase I (before water is exhausted)
if mW > 0
    p = (vAi / vA)^1.4 * pAi;
    Ve = sqrt(2 * (p - Patm) / rhow);
    mdot = Cdd * rhow * At * Ve;
    F = 2 * Cdd * (p - Patm) * At;
    dvdt = Cdd * At * sqrt((2 / rhow) * (pAi* (vAi / vA)^1.4 - Patm));
    dmdt = -mdot;
    dmwdt = -mdot;
    dmadt = 0;
end

% Thermodynamic/mass equations for phase II (after water is exhausted)

pend = (vAi / vB)^1.4 * pAi;
p = pend * (mA / mAi)^1.4;
if mW == 0 && p > Patm
    rhoa = mA / vB;
    T = p / (rhoa * R);
    pcrit = p * (2/2.4)^(1.4/0.4);
    if pcrit > Patm
        pe = pcrit;
        Te = (2 / 2.4) * T;
        rhoe = pcrit / (R * Te);
        Ve = sqrt(1.4 * R * Te);
    else
        pe = Patm;
        Me = sqrt(((p / Patm)^(0.4/1.4) - 1) * 2 / 0.4);
        Te = T * (1 + 0.4 * Me^2 / 2);
        rhoe = Patm / (R * Te);
        Ve = Me * sqrt(1.4 * R * Te);
    end
    mdot = Cdd * rhoe * At * Ve; %%%% REMEMBER that this used to be Cd not Cdd, might need some fixing%%%%
    F = mdot * Ve + At * (pe - Patm);
    dmadt = -mdot;
    dmdt = -mdot;
    dmwdt = 0;
    dvdt = 0;    
end

if p < Patm
    p = Patm;
end

% Thermodynamic/mass equations for phase III (ballistic)
if mW == 0 && p == Patm
    F = 0;
    dmdt = 0;
    dmadt = 0;
    dmwdt = 0;
    dvdt = 0;
end

% Run differential equation results

dzdt = V * sin(theta);
dxdt = V * cos(theta);
% dmdt is provided per the conditional statements above
dVdt = (F - D - mR * g0 * sin(theta)) / mR;
dthdt = (-g0 * cos(theta)) / (V);
% dvdt is provided per the conditional statements above
% dmwdt is provided per the conditional statements above
% dmadt is provided per the conditional statements above

% Conditions for hitting the ground

if z <= 0
    dzdt = 0;
    dxdt = 0;
end

%Check for and correct anomolous conditions

if V < 0.00001
   dthdt = 0; 
end
if -pi/2 > theta
    dthdt = 0;
end

% Assign variables to their alloted slot in the h matrix

dhdt(1) = dzdt; % Change in vertical position
dhdt(2) = dxdt; % Change in horizontal position
dhdt(3) = dmdt; % Change in mass of the entire rocket
dhdt(4) = dVdt; % Change in velocity
dhdt(5) = dthdt; % Change in tragectory angle
dhdt(6) = dvdt; % Change in volume of the air in the bottle
dhdt(7) = dmwdt; % Change in mass of water
dhdt(8) = dmadt; % Change in mass of air

end