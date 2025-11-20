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
SOC_initial = 0.5;
% C-rate for charging/discharging (e.g., 1C = full charge/discharge in 1 hour)
C_rate = NaN;
thermal_mode = NaN;
% Operating temperature in degrees Celsius
temperature = 25;
% Analysis mode: "full-cell" or "half-cell"
cell_analysis_mode = "full-cell";
cycle_startmode = "discharge";

%% Run Impedance
% ---------------- User settings ----------------
freq_min = 200e-6;    % 200 µHz
freq_max = 1e3;       % 1 kHz
numFreqPoints = 65;   % Number of frequency points

I_app_amplitude = 0.05;           % Applied current amplitude [A]

samples_per_period = 100;         % Number of time points per period
num_periods_total  = 40;       % total periods per frequency
num_periods_steady = 10;       % use only last 10 for FFT

% ------------------------------------------------
frequencies = logspace(log10(freq_min), log10(freq_max), numFreqPoints);  % descending order
% frequencies = flip(frequencies);  % optional: make it ascending

% ---------------- Output array ----------------
impedance = zeros(length(frequencies), 1);

for k = 1:numFreqPoints
    f = frequencies(k);
    T = 1/f;
    dt = T / samples_per_period;
    
    % ---- time vector for the frequency ----
    N  = samples_per_period * num_periods_total;   % total samples
    t  = (0:N-1).' * dt;

    % ---- sinusoidal current ----
    applied_current = I_app_amplitude * sin(2*pi*f*t);
        
    % ---- run time-domain simulation ----
    timeFinal = t(end) + 10*eps(t(end));  % just beyond last sample
    timeStep = dt;    % Time step size in seconds

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
    applied_current_vectors(:, k) = applied_current; % Cell array to store vectors

    if numel(voltage) ~= N
        % Assuming v_raw is given at constant step dt_sim = timeStep
        t_sim = (0:numel(voltage)-1).' * dt;
        v = interp1(t_sim, voltage, t, 'linear', 'extrap');
    else
        v = voltage;
    end

    % ---- use only last num_periods_steady periods for FFT ----
    M = samples_per_period * num_periods_steady;  % samples in window
    start_idx = N - M + 1;

    I_win = applied_current(start_idx:end);
    V_win = v(start_idx:end);

    % (optional) remove mean to reduce leakage
    I_win = I_win - mean(I_win);
    V_win = V_win - mean(V_win);

    % ---- FFT and fundamental extraction ----
    N_win = numel(I_win);
    I_fft = fft(I_win) / N_win;
    V_fft = fft(V_win) / N_win;

    % frequency resolution for this window
    df  = 1 / (N_win * dt);
    % index of our excitation frequency:
    idx = round(f / df) + 1;   % +1 because FFT index 1 is DC

    % safety check (optional)
    if idx < 1 || idx > N_win
        warning('Frequency %.4g Hz out of FFT range', f);
        impedance(k) = NaN;
        continue;
    end

    % complex impedance at that frequency
    impedance(k) = V_fft(idx) / I_fft(idx);

    % num_points_last_three = 3 * samples_per_period + 1;
    % start_idx = length(t_eval) - num_points_last_three + 1;
    % current_last = applied_current(start_idx:end);
    % voltage_last = voltage(start_idx:end);
    % 
    % % Apply FFT to cell voltage
    % voltage_fft = fft(voltage_last) / length(voltage_last);
    % current_fft = fft(current_last) / length(current_last);
    % 
    % % [~, idx] = max(abs(current_fft));  % peak at fundamental freq
    % idx = round(f / (1/(num_points*dt))) + 1;
    % 
    % % Compute impedance
    % impedance_tmp = voltage_fft(idx) / current_fft(idx);
    % impedance(i) = impedance_tmp;
   
end

%% Plot Nyquist Diagram
figure
plot(real(impedance),-imag(impedance),'. -k','MarkerSize',10);
axis equal
xlabel('$Z_\mathrm{r}(\omega)$ [m$\Omega$]','Interpreter','latex','FontSize',11);
ylabel('$-Z_\mathrm{j}(\omega)$ [m$\Omega$]','Interpreter','latex','FontSize',11); 

