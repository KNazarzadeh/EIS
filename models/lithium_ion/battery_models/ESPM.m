classdef ESPM
    % ESPM class for Single/Extended Single Particle Model battery simulation.
    %
    % This class configures a battery object using the ESPM framework.
    % It supports both full-cell and half-cell configurations and initializes
    % key parameter sets required for the simulation, including geometry,
    % transport, kinetics, thermodynamics, and electrical behavior.
    %
    methods(Static)
        function ESPMinitialParameters( ...
                battery, ...
                applied_current, ...
                working_electrode, ...
                cycleNum, ...
                timeStep, ...
                temperature)

            % Constructor for the SPM class.
            %
            % Inputs:
            %   battery            - Battery object with base parameters
            %   applied_current    - Time-varying applied current profile
            %   cell_analysis_mode - Either 'full-cell' or 'half-cell'
            %   working_electrode  - Electrode(s) to simulate ('positive', 'negative', or both)
            %
            % This constructor updates the battery object to include SPM-specific
            % parameter structures for geometry, transport, concentration, kinetics,
            % thermodynamics, and electrical behavior.
            
            % -------- Geometric parameters --------
            GeometricParameters.initialGeometricParameters(battery, applied_current, working_electrode, cycleNum);

            % -------- Precompute FVM matrices (electrodes) --------
            for electrode = working_electrode                
                battery.preComputedMatrix.electrode.(electrode) = FVMelectrode.compute_electrodeFVMMatirces( ...
                    battery.GeometricParams, battery.TransportParams, timeStep, electrode);
            end
            
            % -------- Concentration parameters --------
            ConcentrationParameters.initialConcentrationParameters(battery, battery.preComputedMatrix.electrode, ...
                working_electrode, cycleNum);

            % -------- Kinetic parameters --------
            KineticParameters.initialKineticParameters(battery, applied_current, working_electrode, cycleNum);

            % -------- Transport parameters --------
            TransportParameters.initialTransportParameters(battery);

            % -------- Electrolyte FVM matrices (fix: do NOT use undefined `electrode`) --------
            % ------------------------ Full Cell Analysis ------------------------
            battery.preComputedMatrix.electrolyte = FVMelectrolyte.compute_electrolyteFVMMatirces_FullCell( ...
                battery.GeometricParams, battery.TransportParams, battery.ThermalParams, ...
                battery.ConcentrationParams.electrolyte.concentration_initial_elyte, temperature, timeStep);

            % -------- Thermodynamic parameters --------
            ThermodynamicParameters.initialThermodynamicParameters(battery, battery.preComputedMatrix.electrolyte, ...
                working_electrode, cycleNum);

            % -------- Electrical parameters --------
            ElectricalParameters.initialElectricalParameters(battery, applied_current, cycleNum);
        end
    end
end