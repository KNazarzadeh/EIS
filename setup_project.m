% setup_project.m
% Configures the MATLAB path for the V2G project

% Get the root directory of the project
root_dir = fileparts(mfilename('fullpath'));

% Add constant file to the MATLAB path
addpath(fullfile(root_dir, 'datapreprocessing'));

% Add Parameter Category files to the MATLAB path
addpath(fullfile(root_dir, 'parameters', 'DefaultParameters'));
addpath(fullfile(root_dir, 'parameters', 'ConstantParameters'));
addpath(fullfile(root_dir, 'parameters', 'ConcentrationParameters'));
addpath(fullfile(root_dir, 'parameters', 'GeometricParameters'));
addpath(fullfile(root_dir, 'parameters', 'ElectricalParameters'));
addpath(fullfile(root_dir, 'parameters', 'KineticParameters'));
addpath(fullfile(root_dir, 'parameters', 'ThermodynamicParameters'));
addpath(fullfile(root_dir, 'parameters', 'TransportParameters'));
addpath(fullfile(root_dir, 'parameters', 'ThermalParameters'));
addpath(fullfile(root_dir, 'parameters', 'UpdateParameters'));

% Add Simulation folder
addpath(fullfile(root_dir, 'simulation'));

% Add project folders to the MATLAB path
addpath(fullfile(root_dir, 'models', 'lithium_ion', 'battery_models'));
addpath(fullfile(root_dir, 'models', 'lithium_ion', 'aging_models'));

addpath(fullfile(root_dir, 'models', 'lithium_ion', 'submodels', 'FVM'));
addpath(fullfile(root_dir, 'models', 'lithium_ion', 'submodels', 'CCCV'));
addpath(fullfile(root_dir, 'models', 'lithium_ion', 'submodels', 'PackLevel'));
% Add common files
addpath(fullfile(root_dir, 'utils'));

addpath(fullfile(root_dir, 'Test'));
addpath(fullfile(root_dir, 'plots'));

% Add Optimization folder
addpath(fullfile(root_dir, 'Optimization'));

% Optional: Confirm setup
disp('Project paths configured successfully.');