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
                applied_current, ...
                cycleNum)
            % INITIALGEOMETRICPARAMETERS  Populate ElectricalParams efficiently
            %
            %   obj = initialElectricalParameters(obj)
            %

            % --------------------------------------------------------------------
            % 1. Initialise calculator & aging flags
            % --------------------------------------------------------------------
            isLaterAging = (cycleNum > 1);
            calculator = ElectricalParametersCalculator();

            ElecCell = battery.ElectricalParams.cell;
            % --------------------------------------------------------------------
            % 2. Cell domain
            % --------------------------------------------------------------------
            % -------- Cell domain --------
            [ElecCell.OCV_min, ...
             ElecCell.OCV_max] = calculator.compute_ocv_range(battery);

            ElecCell.SOC = [];
            ElecCell.SOC_initial = calculator.compute_initial_SOC(battery, cycleNum);
            ElecCell.SOC(1) = ElecCell.SOC_initial;

            % -------- Full-cell voltage ----------------------------------------
                if isLaterAging
                    ElecCell.voltage = [];
                    ElecCell.cell_capacity = [];

                    total_capacityLoss = ElecCell.total_capacityLoss(end);
                    battery_capacity = ElecCell.battery_capacity_cell(cycleNum-1, end);
                    ElecCell.battery_capacity_cell(cycleNum, 1) = ...
                    battery_capacity - total_capacityLoss;

                    ElecCell.total_capacityLoss =  [];
                    ElecCell.total_capacityLoss(1) = total_capacityLoss;

                    capacityLossLAM_neg = battery.ElectricalParams.electrode.negative.capacityLossLAM_neg(end);
                    battery.ElectricalParams.electrode.negative.capacityLossLAM_neg = [];
                    battery.ElectricalParams.electrode.negative.capacityLossLAM_neg(1) = capacityLossLAM_neg;

                    capacityLossLAM_pos = battery.ElectricalParams.electrode.positive.capacityLossLAM_pos(end);
                    battery.ElectricalParams.electrode.positive.capacityLossLAM_pos = [];
                    battery.ElectricalParams.electrode.positive.capacityLossLAM_pos(1) = capacityLossLAM_pos;
    
                    capacityLossSEI_neg = battery.ElectricalParams.electrode.negative.capacityLossSEI_neg(end);
                    battery.ElectricalParams.electrode.negative.capacityLossSEI_neg = [];
                    battery.ElectricalParams.electrode.negative.capacityLossSEI_neg(1) = capacityLossSEI_neg;

                    sei_film_resistance_neg = battery.ElectricalParams.electrode.negative.sei_film_resistance_neg(end);
                    battery.ElectricalParams.electrode.negative.sei_film_resistance_neg = [];
                    battery.ElectricalParams.electrode.negative.sei_film_resistance_neg(1) = sei_film_resistance_neg;
                else
                    ElecCell.battery_capacity_cell(1, 1) = ...
                        calculator.compute_battery_capacity(battery, cycleNum);

                    battery.ElectricalParams.electrode.negative.sei_film_resistance_neg(1) = ...
                        battery.ElectricalParams.electrode.negative.sei_film_resistance_neg;
                    battery.ElectricalParams.electrode.negative.capacityLossLAM_neg(1) = 0;
                    battery.ElectricalParams.electrode.positive.capacityLossLAM_pos(1) = 0;
                    battery.ElectricalParams.electrode.negative.capacityLossSEI_neg(1) = 0;
                end

                ElecCell.voltage(1) = calculator.compute_cell_voltage(battery, applied_current, 0);
                
                battery.ElectricalParams.cell = ElecCell;
            
            % --------------------------------------------------------------------
            % 4. Model-specific geometry (executed once, after electrodes)
            % --------------------------------------------------------------------
            % Extended Single Particle Model (ESPM) uses basic OCP parameters
        end
    end
end