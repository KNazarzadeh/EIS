clear all;
close all;
clc;
% profile clear                              % reset old data
% profile on -timer real                     % wall-clock timing (good default)
%% Set Input file
databasePath = "C:\Users\k.nazarzadeh\Projects\EIS";
inputFileName = "knazarzadeh_ESPM_Pouch_cells_EIS_2.xlsx";

%% Battery Parameters
% Read Battery Parameters from Excel file
batteryStruct = readExcel(fullfile(databasePath, inputFileName));
battery = Battery(batteryStruct);   % <-- now it's a handle
%% Simulation Parameters
model = "ESPM";

% Initial State of Charge (SOC) as a fraction (0 = 0%, 1 = 100%)
% Example: 0 or "0%" or '0%'
SOC_initial = 0.5;
% C-rate for charging/discharging (e.g., 1C = full charge/discharge in 1 hour)
C_rate = NaN;
thermal_mode = NaN;
% Operating temperature in degrees Celsius
temperature = 25;
% Analysis mode: "full-cell" or "half-cell"
cell_analysis_mode = "full-cell";
cycle_startmode = "charge";

%% Run Impedance
% ---------------- User settings ----------------
freq_min = 200e-6;    % 200 µHz
freq_max = 1e3;       % 1 kHz
numFreqPoints = 65;   % Number of frequency points

I_app_amplitude = 0.05;           % Applied current amplitude [A]

samples_per_period = 100;         % Number of time points per period
num_periods= 10;                  % use only last 10 for FFT

% ------------------------------------------------
frequencies = logspace(log10(freq_min), log10(freq_max), numFreqPoints);  % descending order
% frequencies = flip(frequencies);  % optional: make it ascending

% ---------------- Output array ----------------
impedance = zeros(numFreqPoints, 1);

voltage_vectors = zeros(1000, numFreqPoints);
overpotential_negative_vectors = zeros(1000, numFreqPoints);
overpotential_positive_vectors = zeros(1000, numFreqPoints);
potential_negative_vectors = zeros(1000, numFreqPoints);
potential_positive_vectors = zeros(1000, numFreqPoints);

for k = 1:numFreqPoints
    f = frequencies(k);
    period = 1/f;
    timeStep = period / samples_per_period;  % Time step size in seconds
    
    % ---- time vector for the frequency ----
    N_sample  = samples_per_period * num_periods;   % total samples
    t_eval = (0:(N_sample) - 1) * timeStep;

    % ---- sinusoidal current ----
    applied_current = I_app_amplitude * sin(2*pi*f*t_eval).';
        
    % ---- run time-domain simulation ----
    timeFinal = t_eval(end) + 10*eps(t_eval(end));  % just beyond last sample

    % ---- Run Simulation ----
    simulation(model, ...
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

    % --- Get simulated voltage and align vectors ---
    voltage = battery.ElectricalParams.cell.voltage;
    voltage_vectors(:, k) = voltage;

    overpotential_negative = battery.ThermodynamicParams.electrode.negative.overpotential_neg;
    overpotential_positive = battery.ThermodynamicParams.electrode.positive.overpotential_pos;
    potential_negative = battery.ThermodynamicParams.electrode.negative.potential_neg;
    potential_positive = battery.ThermodynamicParams.electrode.positive.potential_pos;

    overpotential_negative_vectors(:, k) = overpotential_negative;
    overpotential_positive_vectors(:, k) = overpotential_positive;
    potential_negative_vectors(:, k) = potential_negative;
    potential_positive_vectors(:, k) = potential_positive;


    num_points_last_three = 3 * samples_per_period + 1;
    start_idx = length(t_eval) - num_points_last_three + 1;
    current_last = applied_current(start_idx:end);
    voltage_last = voltage(start_idx:end);

    overpotential_negative_last = overpotential_negative(start_idx:end);
    overpotential_positive_last = overpotential_positive(start_idx:end);
    potential_negative_last = potential_negative(start_idx:end);
    potential_positive_last = potential_positive(start_idx:end);

    % Apply FFT to cell voltage
    voltage_fft = fft(voltage_last) / length(voltage_last);
    current_fft = fft(current_last) / length(current_last);

    % Apply FFT to overpotential Negative
    overpotential_negative_fft = fft(overpotential_negative_last) / length(overpotential_negative_last);
    % Apply FFT to overpotential Psitive
    overpotential_positive_fft = fft(overpotential_positive_last) / length(overpotential_positive_last);
    % Apply FFT to potential Negative
    potential_negative_fft = fft(potential_negative_last) / length(potential_negative_last);
    % Apply FFT to potential Positive
    potential_positive_fft = fft(potential_positive_last) / length(potential_positive_last);


    % Get first harmonic (frequency index ≈ 5th in this case)
    [~, idx] = max(abs(current_fft));  % peak at fundamental freq

    % Compute impedance
    impedance_tmp = voltage_fft(idx) / current_fft(idx);
    
    impedance(k) = impedance_tmp;
    impedance_negative_overpotential(k) = overpotential_negative_fft(idx) / current_fft(idx);
    impedance_positive_overpotential(k) = overpotential_positive_fft(idx) / current_fft(idx);
    impedance_negative_potential(k) = potential_negative_fft(idx) / current_fft(idx);
    impedance_positive_potential(k) = potential_positive_fft(idx) / current_fft(idx);

    fprintf('i = %.4f s\n', k);   
end

figure
plot(1000*real(impedance), -1000*imag(impedance), '-ok' , 'LineWidth' , 2 );
% axis equal
hold on

figure
plot(1000*real(impedance_negative_overpotential), -1000*imag(impedance_negative_overpotential), '-ok' , 'LineWidth' , 2 );
% axis equal
hold on

figure
plot(1000*real(impedance_positive_overpotential), -1000*imag(impedance_positive_overpotential), '-ok' , 'LineWidth' , 2 );
% axis equal
hold on


fprintf('Total simulation time: %.2f seconds (%.2f minutes)\n', battery.SimulationParams.('elapsedTime'), battery.SimulationParams.('elapsedTime')/60);

disp("------------------------")
save("batteryStruct_2.mat", "batteryStruct");