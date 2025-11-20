classdef KineticParameters
    % Extract and organize kinetic parameters for battery models
    %
    % Description:
    % -----------
    % This class extracts kinetic parameters from a specified battery type
    % (e.g., NMC811_chen2020) and organizes them into domains (cell, electrode,
    % electrolyte, separator, SEI) based on the electrochemical model type
    % (e.g., SPM, ESPM, DFN).
    %
    
    methods(Static)
        % ---------------------------------------------------------------
        function initialKineticParameters( ...
                battery, ...
                applied_current, ...
                working_electrode)
            % Initial and organize kinetic parameters by domain
            %
            % Description:
            % ------------
            % extracts kinetic parameters, organizing them into substructures:
            % cell, electrode (negative, positive, including particle parameters),
            % electrolyte, separator, and SEI based on the battery type and
            % the specified electrochemical model (SPM, ESPM, etc.).
            %
            % Outputs:
            % --------
            % kineticParams: Updated object with kinetic_params structure populated

            % --------------------------------------------------------------------
            % 1. Initialise calculator
            % --------------------------------------------------------------------
            calculator = KineticParametersCalculator();

            % --------------------------------------------------------------------
            % 2. Main electrode loop
            % --------------------------------------------------------------------
            % -------- Electrode domain --------
            for electrode = working_electrode
                % ---- Prefix for electrodes ('neg', 'pos') ----------------------
                prefix = electrode{1}(1:3);
                
                % ---- local cache
                Kinet = battery.KineticParams.electrode.(electrode);
                CancE = battery.ConcentrationParams.electrode.(electrode);
                Thrmody = battery.ThermodynamicParams.electrode.(electrode);
                CancElyt = battery.ConcentrationParams.electrolyte;

                % ---- Molar Ionic Flux ----------------------------------------
                Thrmody.(['potential_timederivation_' prefix])(1, :) = 0;
                Kinet.(['molar_ionic_flux_' prefix])(1) = ...
                    calculator.compute_molar_ionic_flux(battery, applied_current(1), ...
                    Thrmody.(['potential_timederivation_' prefix])(1, :), ...
                    electrode);

                % ---- Reaction Rate Coefficient --------------------------------
                Kinet.(['reaction_rate_coefficient_' prefix]) = ...
                    calculator.compute_reaction_rate_coefficient(battery, electrode);

                % ---- Exchange Current Density ----------------------------------------
                battery.KineticParams.electrode.(electrode) = Kinet;
                Kinet = battery.KineticParams.electrode.(electrode);
                Kinet.(['exchange_current_density_' prefix])(1) = ...
                    calculator.compute_exchange_current_density(battery, ...
                    CancE.(['surface_concentration_' prefix])(1), ...
                    CancElyt.(['concentration_mean_elyte_' prefix])(1), ...
                    electrode);

                battery.KineticParams.electrode.(electrode) = Kinet;
            end
            
            % --------------------------------------------------------------------
            % 4. Model-specific Kinetic (executed once, after electrodes)
            % --------------------------------------------------------------------
            % ESPM Model:Extended Single Particle Model (ESPM)
            battery.KineticParams.separator.molar_ionic_flux_sep = 0;
        end
    end
end