function simulation(model, ...
    battery, ...
    SOC_initial, ...
    applied_current, ...
    C_rate, ...
    cycle_startmode, ...
    temperature, ...
    thermal_mode, ...
    cell_analysis_mode, ...
    timeFinal, ...
    timeStep)

    %% Validate battery structure
    if isempty(battery), error('Battery parameter structure is required.'); end

    %% Set default values for optional inputs
    if isnan(SOC_initial), SOC_initial = 1; end
    if isnan(timeFinal), timeFinal = 3600; end
    if isnan(timeStep), timeStep = 1; end
    if isnan(temperature) || isempty(temperature)
        temperature = 298.15;   % Default to 25°C in Kelvin
    else
        temperature = temperature + 273.15; % Convert °C to Kelvin
    end
    %% Validate SOC_initial configuration
    % Check if value is a string and contains '%'
    if ischar(SOC_initial) || isstring(SOC_initial)
        if contains(SOC_initial, '%')
            valueStr = strrep(SOC_initial, '%', ''); % Remove %
            SOC_initial = str2double(valueStr) / 100;    % Convert to number
        else
            SOC_initial = str2double(SOC_initial);          % Convert directly
        end
    end

    %% working electrode configuration
    if strcmp(cell_analysis_mode, 'full-cell') || strcmp(cell_analysis_mode, 'full')
        cell_analysis_mode = 'full-cell';
        if ~exist('working_electrode','var') || isempty(working_electrode) || (isnumeric(working_electrode) && any(isnan(working_electrode)))
            working_electrode = ["negative", "positive"];
        end

    elseif strcmp(cell_analysis_mode, 'half-cell') || strcmp(cell_analysis_mode, 'half')
        cell_analysis_mode = 'half-cell';
        switch lower(working_electrode)
            case {'anode', 'negative'}
                % Run half-cell simulation with anode
                working_electrode = 'negative';
            case {'cathode', 'positive'}
                % Run half-cell simulation with cathode
                working_electrode = 'positive';
            otherwise
                error('working_electrode must be ''anode'', ''negative'', ''cathode'', or ''positive''.');
        end
    end    
   %% Store simulation parameters in battery struct

    % Simulation start time in seconds
    timeStart = 1;
    battery.SimulationParams = struct( ...
        'SOC_initial', SOC_initial, ...
        'C_rate', C_rate, ...
        'timeStart', timeStart, ...
        'timeFinal', timeFinal, ...
        'timeStep', timeStep, ...
        'cycle_startmode', cycle_startmode, ...
        'temperature', temperature, ...
        'thermal_mode', thermal_mode, ...
        'cell_analysis_mode', cell_analysis_mode);
            
    timeSize = round(timeFinal/timeStep);    
    battery.SimulationParams.timeSize = timeSize;

    %% Run simulation for selected model
    % Determine applied current using C-rate if needed
    applied_current = DataPreprocessing.preprocessCurrent(battery, applied_current, C_rate, cycle_startmode);

    switch model
        case "ESPM"
            ESPM.ESPMinitialParameters(battery, applied_current, working_electrode, timeStep, temperature);
            ESPMsimulation.runESPMsimulation(battery, battery.preComputedMatrix, applied_current, timeSize, timeStep);
    end

end