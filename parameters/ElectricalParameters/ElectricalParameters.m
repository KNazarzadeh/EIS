classdef ElectricalParameters
    % Extract and organize thermodynamic parameters (OCP) for battery models
    %
    % Description:
    % -----------
    % This class extracts thermodynamic parameters, specifically Open-Circuit Potential (OCP)
    % from a specified battery type (e.g., NMC811_chen2020) and organizes them into domains (cell, electrode,
    % electrolyte, separator, sei) based on the electrochemical model type
    % (e.g., SPM, ESPM, DFN).
    %
    
    methods(Static)
        % ---------------------------------------------------------------
        function initialElectricalParameters( ...
                battery, ...
                applied_current)
            % INITIALGEOMETRICPARAMETERS  Populate ElectricalParams efficiently
            %
            %   obj = initialElectricalParameters(obj)
            %
            % --------------------------------------------------------------------
            % 1. Initialise calculator
            % --------------------------------------------------------------------
            calculator = ElectricalParametersCalculator();

            ElecCell = battery.ElectricalParams.cell;
            % --------------------------------------------------------------------
            % 2. Cell domain
            % --------------------------------------------------------------------
            % -------- Cell domain --------
            [ElecCell.OCV_min, ...
             ElecCell.OCV_max] = calculator.compute_ocv_range(battery);

            ElecCell.SOC_initial = calculator.compute_initial_SOC(battery);
            % -------- Full-cell voltage ----------------------------------------
            ElecCell.battery_capacity_cell(1, 1) = ...
                calculator.compute_battery_capacity(battery);

            ElecCell.voltage(1) = calculator.compute_cell_voltage(battery, applied_current, 0);
                
                battery.ElectricalParams.cell = ElecCell;
            
            % --------------------------------------------------------------------
            % 4. Model-specific geometry (executed once, after electrodes)
            % --------------------------------------------------------------------
            % Extended Single Particle Model (ESPM) uses basic OCP parameters
        end
    end
end