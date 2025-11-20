%% Optimization script

clear; clc; close all;
warning('off','all');

% -------------------------------------------------------------------------
% 1) Load battery and measurement data
% -------------------------------------------------------------------------
load("battery.mat", "battery");   % battery is a handle object

% Each Data_mea* is [time, voltage]
Data_mea1 = [MeasurementsREUSE.Low(1:4730, 2),      MeasurementsREUSE.Low(1:4730, 3)];
Data_mea2 = [MeasurementsREUSE.Med(1:1307, 2),      MeasurementsREUSE.Med(1:1307, 3)];
Data_mea4 = [MeasurementsREUSE.MedHigh(:, 2),       MeasurementsREUSE.MedHigh(:, 3)];

% Simulation / pack params (set these to what you already use)
dischargeC_rate = 0.5;       % <-- or your actual value
temperature     = 25;        % [°C]
cycleTotalNum   = 1;         % or more, if you fit over multiple cycles
parallelCellNum = 3;
seriesCellNum   = 96;

% -------------------------------------------------------------------------
% 2) Freeze the *initial* battery as a struct snapshot
% -------------------------------------------------------------------------
batteryStruct0 = battery.toStruct();   % VERY IMPORTANT

% -------------------------------------------------------------------------
% 3) Bounds and initial guess for 6 parameters
% -------------------------------------------------------------------------
lb = [1e-15,  1e-15,  28000,  48000,  4e-12,  4e-12];
ub = [3e-13,  4e-13,  34000,  52000,  1e-9,   1e-9];

x0 = [ ...
    battery.TransportParams.electrode.negative.diffusion_coefficient_reference_neg; ...
    battery.TransportParams.electrode.positive.diffusion_coefficient_reference_pos; ...
    battery.ConcentrationParams.electrode.negative.concentration_maximum_neg; ...
    battery.ConcentrationParams.electrode.positive.concentration_maximum_pos; ...
    battery.KineticParams.electrode.negative.reaction_rate_coefficient_reference_neg; ...
    battery.KineticParams.electrode.positive.reaction_rate_coefficient_reference_pos];

% -------------------------------------------------------------------------
% 4) Objective over multiple operating conditions
% -------------------------------------------------------------------------
objective = @(x) ...
    VoltageError_ESPM(x, batteryStruct0, Data_mea1, dischargeC_rate, temperature, cycleTotalNum, parallelCellNum, seriesCellNum) + ...
    VoltageError_ESPM(x, batteryStruct0, Data_mea2, dischargeC_rate, temperature, cycleTotalNum, parallelCellNum, seriesCellNum) + ...
    VoltageError_ESPM(x, batteryStruct0, Data_mea4, dischargeC_rate, temperature, cycleTotalNum, parallelCellNum, seriesCellNum);

% -------------------------------------------------------------------------
% 5) Optimize
% -------------------------------------------------------------------------
opts = optimset('fminsearch');
opts.Display      = 'iter';
opts.MaxIter      = 100;
opts.MaxFunEvals  = 5000;
opts.TolX         = 1e-6;
opts.TolFun       = 1e-6;

objective_wrapped = @(x) objective_with_bounds(x, lb, ub, ...
    batteryStruct0, Data_mea1, Data_mea2, Data_mea4, ...
    dischargeC_rate, temperature, cycleTotalNum, parallelCellNum, seriesCellNum);

[X, fval, exitflag, output] = fminsearch(objective_wrapped, x0, opts);

% -------------------------------------------------------------------------
% 6) Apply optimized parameters to a battery for later use
% -------------------------------------------------------------------------
battery_opt = Battery(batteryStruct0);

battery_opt.TransportParams.electrode.negative.diffusion_coefficient_reference_neg = X(1);
battery_opt.TransportParams.electrode.positive.diffusion_coefficient_reference_pos = X(2);
battery_opt.ConcentrationParams.electrode.negative.concentration_maximum_neg       = X(3);
battery_opt.ConcentrationParams.electrode.positive.concentration_maximum_pos       = X(4);
battery_opt.KineticParams.electrode.negative.reaction_rate_coefficient_reference_neg = X(5);
battery_opt.KineticParams.electrode.positive.reaction_rate_coefficient_reference_pos = X(6);

% You can save this:
% save('battery_optimized.mat','battery_opt');

%% ------------------------------------------------------------------------
% Local function(s) – must come AFTER the main script in a .m file
% -------------------------------------------------------------------------
function err = objective_with_bounds(x, lb, ub, ...
    batteryStruct0, Data_mea1, Data_mea2, Data_mea4, ...
    dischargeC_rate, temperature, cycleTotalNum, parallelCellNum, seriesCellNum)

    % Clamp x into [lb, ub]
    x_clamped = max(min(x, ub), lb);

    % Sum error over all operating conditions
    err = ...
        VoltageError_ESPM(x_clamped, batteryStruct0, Data_mea1, dischargeC_rate, temperature, cycleTotalNum, parallelCellNum, seriesCellNum) + ...
        VoltageError_ESPM(x_clamped, batteryStruct0, Data_mea2, dischargeC_rate, temperature, cycleTotalNum, parallelCellNum, seriesCellNum) + ...
        VoltageError_ESPM(x_clamped, batteryStruct0, Data_mea4, dischargeC_rate, temperature, cycleTotalNum, parallelCellNum, seriesCellNum);
end

function err = VoltageError_ESPM(x, batteryStruct0, data_mea, ...
                                 dischargeC_rate, temperature, ...
                                 cycleTotalNum, parallelCellNum, seriesCellNum)

    % === Rebuild a fresh battery object from *initial* struct ===
    battery = Battery(batteryStruct0);

    % === Update the 6 parameters from optimization vector x ===
    battery.TransportParams.electrode.negative.diffusion_coefficient_reference_neg = x(1);
    battery.TransportParams.electrode.positive.diffusion_coefficient_reference_pos = x(2);

    battery.ConcentrationParams.electrode.negative.concentration_maximum_neg = x(3);
    battery.ConcentrationParams.electrode.positive.concentration_maximum_pos = x(4);

    battery.KineticParams.electrode.negative.reaction_rate_coefficient_reference_neg = x(5);
    battery.KineticParams.electrode.positive.reaction_rate_coefficient_reference_pos = x(6);

    try
        % === Run your existing ESPM simulation ===
        SOC_initial       = 1;
        applied_current   = NaN;          % use C-rate based profile
        current_mode      = "CCCV";
        cycle_startmode   = "discharge";
        aging_mode        = 1;
        timeFinal         = 2e5;
        timeStep          = 1;
        dischargeTimeStep = 1;
        chargeTimeStep    = 1;
        CVTimeStep        = 1;
        restDuration      = 300;
        CCCVstatus        = 1;
        current_cutoff    = 0.05;

        simulation( ...
            battery, ...
            SOC_initial, ...
            applied_current, ...
            current_mode, ...
            dischargeC_rate, ...
            dischargeC_rate, ...
            cycle_startmode, ...
            temperature, ...
            aging_mode, ...
            timeFinal, ...
            timeStep, ...
            dischargeTimeStep, ...
            chargeTimeStep, ...
            CVTimeStep, ...
            restDuration, ...
            CCCVstatus, ...
            current_cutoff, ...
            cycleTotalNum, ...
            parallelCellNum, ...
            seriesCellNum);

        % === Extract simulated time and voltage ===
        sim_t = battery.SimulationParams.timeVector(:);
        sim_V = battery.ElectricalParams.cell.voltage(:);

        % === Interpolate to measurement times ===
        t_mea = data_mea(:,1);
        V_mea = data_mea(:,2);

        V_interp = interp1(sim_t, sim_V, t_mea, 'linear', 'extrap');

        % === RMS error in mV ===
        err = rms(V_interp - V_mea) * 1000;

    catch ME
        warning("Simulation failed in VoltageError_ESPM: %s", ME.message);
        err = 1e6;   % heavy penalty for failures
    end
end

disp("------------------------")
disp("Completed!")





