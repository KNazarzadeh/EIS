clear all;
close all;
clc;
% profile clear                              % reset old data
% profile on -timer real                     % wall-clock timing (good default)
%% Set Input file
databasePath = "C:\Users\k.nazarzadeh\Projects\EIS";

inputFileName = "knazarzadeh_ESPM_Pouch_cells_EIS.xlsx";
%% Battery Parameters
% Read Battery Parameters from Excel file
batteryStruct = readExcel(fullfile(databasePath, inputFileName));
battery = Battery(batteryStruct);   % <-- now it's a handle
%% Simulation Parameters
model = "ESPM";

% Initial State of Charge (SOC) as a fraction (0 = 0%, 1 = 100%)
% Example: 0 or "0%" or '0%'
SOC_initial = 1;

% C-rate for charging/discharging (e.g., 1C = full charge/discharge in 1 hour)
C_rate = NaN;

thermal_mode = NaN;

% Constant applied current in Amperes; NaN if using C-rate or experimental profile
applied_current = 0.01;

% Operating temperature in degrees Celsius
temperature = 25;

% Simulation end time in seconds
timeFinal = 2e5;
% Time step size in seconds
timeStep = 1;
% Analysis mode: "full-cell" or "half-cell"
cell_analysis_mode = "full-cell";
cycle_startmode = "discharge";

%% Run Simulation
simulation( ...
    model, ...
    battery, ...
    SOC_initial, ...
    applied_current, ...
    C_rate, ...
    cycle_startmode, ...
    temperature, ...
    thermal_mode, ...
    cell_analysis_mode, ...
    timeFinal, ...
    timeStep);


