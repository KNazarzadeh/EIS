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
applied_current = NaN;

% Operating temperature in degrees Celsius
temperature = 25;

% Simulation end time in seconds
timeFinal = 2e5;
% Time step size in seconds
timeStep = 1;
% Analysis mode: "full-cell" or "half-cell"
cell_analysis_mode = "full-cell";
cycle_startmode = "discharge";

%% Run Impedance
I_app_amplitude = 0.01;           % Applied current amplitude [A]

% freq_min = 200e-6;    % 200 µHz
% freq_max = 1e3;       % 1000 Hz (1 kHz)

numFreqPoints = 15;            % Number of frequency points

% frequencies = logspace(log10(freq_min), log10(freq_max), numFreqPoints);  % descending order

frequencies = logspace(log10(0.1), log10(10000), numFreqPoints);  % descending order

samples_per_period = 400;         % Number of time points per period
num_periods = 13;   

impedance = zeros(length(frequencies), 1);

for i = 1:length(frequencies)

    f0  = frequencies(i);            % fundamental test frequency [Hz]
    P   = num_periods;                % number of recorded periods
    period = 1/f0;                       % period [s]
    Ts  = period / samples_per_period;  % sampling period [s]
    fs  = 1/Ts;                       % sampling frequency [Hz]
    T   = P * period;                   % total record length [s]
    t   = (0:Ts:(T - Ts));          % time vector (N x 1), no duplicate endpoint
    N   = numel(t);                   % number of samples

    % --- Excitation current (single sine) ---
    applied_current = I_app_amplitude * sin(2*pi*f0*t)';

    % --- Run electrochemical model (ensure time settings are passed) ---
    timeStep = Ts;       % <- give your simulator the same Ts
    timeFinal = T;       % <- simulate exactly P periods

    % Call your simulation model:
    % Run Simulation
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
    timeVector = battery.SimulationParams.timeVector;  % Time axis [s]

    % --- Keep steady-state only: use an EXACT integer number of periods ---
    K = 4;                                 % keep last K periods (3–5 is fine)
    Ns = samples_per_period;
    Nkeep = K*Ns;                          % *** no +1 ***
    start_idx = length(t) - Nkeep + 1;
    idx = start_idx:length(t);

    t_keep = t(idx).';                     % column
    i_keep = applied_current(idx);
    v_keep = voltage(idx);

    % --- Lock-in (least-squares) extraction of the fundamental ---
    % model: y ≈ a*sin(ωt) + b*cos(ωt) + c
    X = [sin(2*pi*f0*t_keep) cos(2*pi*f0*t_keep) ones(numel(t_keep),1)];
    theta_i = X \ i_keep;                  % [a_i; b_i; c_i]
    theta_v = X \ v_keep;                  % [a_v; b_v; c_v]

    % complex phasors (Re = cos term, Im = -sin term)
    I = theta_i(2) - 1j*theta_i(1);
    V = theta_v(2) - 1j*theta_v(1);

    impedance(i) = V / I;
   
end

%% Plot Nyquist Diagram
figure
plot(real(impedance),-imag(impedance),'. -k','MarkerSize',10);
axis equal
xlabel('$Z_\mathrm{r}(\omega)$ [m$\Omega$]','Interpreter','latex','FontSize',11);
ylabel('$-Z_\mathrm{j}(\omega)$ [m$\Omega$]','Interpreter','latex','FontSize',11); 

