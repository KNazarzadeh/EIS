function outputFilePath = buildOutputPath(battery, simulation_id)
% buildOutputPath generates a simulation output file name based on battery parameters.
%
% Inputs:
%   battery         - struct containing SimulationParams, including:
%                     .model (string), .cell_analysis_mode (string),
%                     .experiment (NaN or number), .C_rate (string or number)
%   simulation_id   - optional ID string or number to distinguish runs
%
% Output:
%   outputFileName  - constructed file name as a string

    % Handle optional simulation_id
    if nargin < 2 || isempty(simulation_id)
        simulation_id = "1";  % Default
    else
        simulation_id = string(simulation_id);
    end
    
    % Extract fields
    model = string(battery.SimulationParams.model);
    cell_analysis_mode = string(battery.SimulationParams.cell_analysis_mode);
    experiment = battery.SimulationParams.experiment;

    % Current date string
    dateStr = string(datetime('now','Format','yyyyMMdd'));

    % Determine simulation type and C-rate ID
    if isnan(experiment)
        simulation_data = "C-rate";
        if isfield(battery.SimulationParams, 'C_rate')
            cRate_id = string(battery.SimulationParams.C_rate);
        else
            cRate_id = "";
        end
    else
        simulation_data = "experiment";
        cRate_id = "";  % Not applicable
    end

    if isempty(simulation_id)
        simulation_id = "1";  % Default ID if not provided
    else
        simulation_id = string(simulation_id);
    end

    % Build file name based on available parts
    parts = ["simulation", model, cell_analysis_mode, simulation_data];

    if cRate_id ~= ""
        parts(end+1) = cRate_id;
    end

    parts(end+1) = dateStr;
    parts(end+1) = simulation_id;

    % Join parts with underscores
    outputFileName = strjoin(parts, "_");


    % Build directory path
    baseDir = pwd;  % Current working directory
    outputDir = fullfile(baseDir, "output", model, cell_analysis_mode, dateStr);

    % Create folder structure if it doesn't exist
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end

    outputFilePath = fullfile(outputDir, outputFileName);
end