function [battery, applied_current_next, simulationEnd] = CC_CV(battery, time, applied_current, Voltage, parallelCellNum)
% CC–CV controller (amps-only).
% Holds V at OCV_max in CV using a PI controller.
%
% Expects in battery.SimulationParams:
%   dischargeCurrent   [A]   (positive magnitude)
%   chargeCurrent      [A]   (positive magnitude)
%   applied_current_min[A]   (CV termination threshold, positive)
%   dischargeTimeStep  [s]
%   chargeTimeStep     [s]
%   CVTimeStep         [s]
% Creates/uses:
%   CCCVstatus         0=CC discharge, 1=CC charge, 2=CV
%   dischargeTimeEnd   [index]
%   chargeTimeEnd      [index]
%   error             previous voltage error
%   errorIntegral     integral of error


    % -------- PI (Proportional–Integral) Controller gains (CV stage) ----
    % used during CV to hold the terminal voltage at maximum Voltage.
    % Proportional gain (Pg)
    % Scales the controller’s reaction directly to the instantaneous error.
    % If voltage is below the setpoint maximum Voltage, a higher Pg makes the controller immediately 
    % increase charging current more aggressively.
    % Too high → oscillations or overshoot.
    % Too low → sluggish response.
    % scales the immediate reaction to error changes (units ≈ A/V).
    Proportional_gain = 5; 

    % Integral gain (Ig):
    % Scales the reaction to the accumulated error over time.
    % Useful to remove steady-state offsets (ensures the voltage really sits at maximum Voltage, not slightly below it).
    % Too high → controller “winds up” and oscillates.
    % Too low → controller takes a long time to eliminate small errors.
    % scales the accumulated reaction over time (units ≈ A/(V·s)).
    Integral_gain = 50; 

    % -------- Set simulation end ----
    % flag to indicate whether the simulation should end after the iteration.
    % Has to be set at every iteration
    simulationEnd = 0; 
    current_cutoff = battery.AgingParams.('current_cutoff');

    % --------- Stage transitions ---------
    Voltage_min = battery.ElectricalParams.cell.('lower_voltage_cell');
    Voltage_max = battery.ElectricalParams.cell.('upper_voltage_cell');
    
    nominal_capacity = battery.ElectricalParams.cell.('nominal_capacity_cell');

    % --------------- Defaults for rest ---------------
    if ~isfield(battery.SimulationParams,'restDuration_s')
                battery.SimulationParams.('restDuration_s') = 300; % 5 minutes
    end
    if ~isfield(battery.SimulationParams,'restTimeStep')
        battery.SimulationParams.('restTimeStep') = 1;     % 1-second steps during rest
    end

    % --------------- Initialize  ---------------
    if time == 1
        % Determine start direction
        % 0: CC discharge, 1: CC charge, 2: CV
        if isfield(battery.SimulationParams,'cycle_startmode')
            if strcmpi(battery.SimulationParams.('cycle_startmode'), "charge")
                battery.AgingParams.CCCVstatus = 1; 
                startmodeIsCharge = true;
    
            elseif strcmpi(battery.SimulationParams.('cycle_startmode'), "discharge")
                battery.AgingParams.CCCVstatus = 0; 
                startmodeIsCharge = false;
    
            else % fallback to applied_current sign
                startmodeIsCharge = (applied_current >= 0);
                battery.AgingParams.CCCVstatus = startmodeIsCharge * 1 + (~startmodeIsCharge) * 0;
            end
        else
            startmodeIsCharge = (applied_current >= 0);
            battery.AgingParams.CCCVstatus = startmodeIsCharge * 1 + (~startmodeIsCharge) * 0;
        end

        battery.SimulationParams.('dischargeTimeEnd') = 1e6;
        battery.SimulationParams.('chargeTimeEnd') = 1e6;
        % PI memory
        battery.AgingParams.('Ierror') = 0;          % previous voltage error (for PI)
        battery.AgingParams.('Perror') = 0;  % integral error
        % Rest bookkeeping
        battery.SimulationParams.('t_rest_start') = -1;
        battery.SimulationParams.('t_rest_end')   = -1;
        
        % Store what we started with
        battery.SimulationParams.startedWithCharge = startmodeIsCharge;

    end

    % ---------------- State transitions ----------------
    % 0: CC discharge, 1: CC charge, 2: CV
    switch battery.AgingParams.('CCCVstatus')
        case 0  % CC discharge
            if (Voltage <= Voltage_min)
                battery.SimulationParams.('dischargeTimeEnd') = time;
                battery.AgingParams.('CCCVstatus') = 3;  % go to rest after discharge
                battery.SimulationParams.('t_rest_start') = time;
                battery.SimulationParams.('t_rest_end')   = time + ...
                    ceil(battery.SimulationParams.('restDuration_s') / battery.SimulationParams.('restTimeStep'));
            end
        case 1
            if Voltage >= Voltage_max
                battery.AgingParams.('CCCVstatus') = 2;                 % switch to CV hold
            end

        case 2  % CV -> Rest-after-charge when taper current reached
            if (applied_current < current_cutoff) && battery.AgingParams.('CCCVstatus') == 2
                battery.SimulationParams.('chargeTimeEnd') = time;
                battery.AgingParams.('CCCVstatus') = 4;  % rest after charge
                battery.SimulationParams.('t_rest_start') = time;
                battery.SimulationParams.('t_rest_end')   = time + ...
                    ceil(battery.SimulationParams.('restDuration_s') / battery.SimulationParams.('restTimeStep'));
            end

            if applied_current == current_cutoff && battery.AgingParams.CCCVstatus == 2
                simulationEnd = 1;
            end

        case 3 % Rest after discharge -> CC charge when rest time elapsed
            if time >= battery.SimulationParams.('t_rest_end')
                battery.AgingParams.('CCCVstatus') = 1;  % CC charge
            end

        case 4  % Rest after charge -> stop when rest time elapsed
            if time >= battery.SimulationParams.('t_rest_end')
                simulationEnd = 1;
            end
    end
    
    % --------- Output current ---------
    switch battery.AgingParams.('CCCVstatus')
        case 0  % CC discharge
            applied_current_next = -battery.SimulationParams.('dischargeC_rate') * nominal_capacity;
            applied_current_next = applied_current_next/parallelCellNum;
            battery.SimulationParams.('timeStep') = battery.SimulationParams.('dischargeTimeStep');
    
        case 1  % CC charge
            applied_current_next = battery.SimulationParams.('chargeC_rate') * nominal_capacity;
            applied_current_next = applied_current_next/parallelCellNum;
            battery.SimulationParams.('timeStep') = battery.SimulationParams.('chargeTimeStep');
    
        case 2  % CV hold with PI on voltage error
    
            % Reduce Voltage to a scalar target
            % Integral error
            battery.AgingParams.('Ierror')(time+1) = Voltage_max - Voltage; 
   
            % Calculating time derivative of the error,
            % dIerror/dt ==> Ierror(t+1) - Ierror(t) / timeStep
            battery.AgingParams.('Perror')(time+1) = ...
                (battery.AgingParams.('Ierror')(time+1) - battery.AgingParams.('Ierror')(time)) / ...
                 battery.SimulationParams.('timeStep');
    
            % This is the velocity form of a PI controller on voltage error, updating current:
            % ΔI=Z. Δt(Ki*IntegralError(t+1)+Kp​*ProportionalError(t+1))​
            % I(t+1) ​= I(t) ​+ dt​ * (Ki * Ierror(t+1)+ Kp​ * Perror(t+1)).
            % If the pack is below the target voltage (error>0) the term KieΔt increases charge current to raise voltage.
            % If the error is changing fast, the KpΔe term damps/accelerates the change to avoid overshoot.
            applied_current_next = applied_current + ...
                battery.SimulationParams.('timeStep') * ...
                (Integral_gain* battery.AgingParams.('Ierror')(time+1) + ...
                Proportional_gain * battery.AgingParams.('Perror')(time+1));
            
            if isfield(battery.SimulationParams,'CVTimeStep')
                battery.SimulationParams.('timeStep') = battery.SimulationParams.('CVTimeStep');
            else
                battery.SimulationParams.('timeStep') = 1;
            end

        case 3  % Rest after discharge
            applied_current_next = 0;
            battery.SimulationParams.('timeStep') = battery.SimulationParams.('restTimeStep');

        case 4  % Rest after charge
            applied_current_next = 0;
            battery.SimulationParams.('timeStep') = battery.SimulationParams.('restTimeStep');
    end
end
