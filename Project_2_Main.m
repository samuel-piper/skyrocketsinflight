function [TestResults, actualDistances, conditionsForDistances, AvarResults, RvarResults, PvarResults, CvarResults] = Project_2_Main()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Purpose:  Provide a one stop shop for execution of all code required for
%               project 2.
%
% How to Call: [TestResults, actualDistances, conditionsForDistances, AvarResults, RvarResults, PvarResults, CvarResults] = Project_2_Main();
%  *******IT TAKES A LITTLE WHILE TO RUN, THERE'S A LOT GOING ON******
%
% Inputs:   None
%
% Outputs:  None
%
% Assumptions:  None
%
% Created: 10/31/16
% Modified: 12/02/16
% Author: a2f341a79180
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

% Run all of the programs

[~, Htest] = Project_2_Test();
[~, ~, actualDistances, conditionsForDistances] = Project_2_Var();
[Avarmax, Rvarmax, Pvarmax, Cvarmax] = Project_2_Sensitivity();
[~, ~] = Project_2_GraphWins();

% Process the results

TestResults = max(Htest(:,1:2));
AvarResults = max(Avarmax)-min(Avarmax);
RvarResults = max(Rvarmax)-min(Rvarmax);
PvarResults = max(Pvarmax)-min(Pvarmax);
CvarResults = max(Cvarmax)-min(Cvarmax);

% Print the results

fprintf('The test rocket landed at %.3g meters and had a maximum height of %.3g meters.\n\n',TestResults(2),TestResults(1))
fprintf('The variance for each parameter (varied from 0.75 to 1.5 times original value) is as follows:\nLaunch Angle - %.3g (m) Height, %.3g (m) Distance\nVolume Ratio - %.3g (m) Height, %.3g (m) Distance\nLaunch Pressure - %.3g (m) Height, %.3g (m) Distance\nDrag Coefficient - %.3g (m) Height, %.3g (m) Distance\n',AvarResults(1),AvarResults(2),RvarResults(1),RvarResults(2),PvarResults(1),PvarResults(2),CvarResults(1),CvarResults(2))
end