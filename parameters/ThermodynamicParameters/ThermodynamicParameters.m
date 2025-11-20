classdef ThermodynamicParameters
    % Extract and organize thermodynamic parameters for battery models
    %
    % Description:
    % -----------
    % This class extracts thermodynamic parameters, from a specified battery type and 
    % organizes them into domains (cell, electrode,
    % electrolyte, separator, sei) based on the electrochemical model type
    % (e.g., SPM, ESPM, DFN).
    %
    
    methods(Static)
        % ---------------------------------------------------------------
        function initialThermodynamicParameters( ...
                battery, ...
                preComputedMatrix, ...
                working_electrode, ...
                cycleNum)
            % Extract and organize thermodynamic parameters by domain
            %
            % Description:
            % ------------
            % Extracts thermodynamic parameters, organizing them into
            % substructures: cell, electrode (negative, positive, including particle parameters),
            % electrolyte, separator, and sei. Parameters are tailored to the specified
            % electrochemical model (SPM, ESPM, DFN).
            %
            % Outputs:
            % --------
            % thermodynamicParams: Updated object with thermodynamicParams structure populated

            % --------------------------------------------------------------------
            % 1. Initialise calculator & aging flags
            % --------------------------------------------------------------------
            isLaterAging = (cycleNum > 1);
            boundaryPotentialElyteDiscreteTime = preComputedMatrix.boundaryPotentialElyteDiscreteTime;
            surfacePotentialElyteMatrix = preComputedMatrix.surfacePotentialElyteMatrix;
            diffusionPotentialElyteMatrix = preComputedMatrix.diffusionPotentialElyteMatrix;

            calculator = ThermodynamicParametersCalculator();

            % --------------------------------------------------------------------
            % 2. Model-specific geometry (executed once, after electrodes)
            % --------------------------------------------------------------------       
            % -------- ESPM: compute electrolyte potential --------
            if isLaterAging
                battery.ThermodynamicParams.electrolyte.potential_elyte = []; 
                battery.ThermodynamicParams.electrolyte.potential_mean_elyte_neg = [];
                battery.ThermodynamicParams.electrolyte.potential_mean_elyte_pos = [];
            end

            battery.ThermodynamicParams.electrolyte.potential_elyte(:, 1) = ...
                calculator.compute_electrolyte_potential( ...
                boundaryPotentialElyteDiscreteTime, ...
                surfacePotentialElyteMatrix, ...
                diffusionPotentialElyteMatrix, ...
                battery.KineticParams.electrode.negative.molar_ionic_flux_neg(1), ...
                battery.KineticParams.electrode.positive.molar_ionic_flux_pos(1), ...
                battery.ConcentrationParams.electrolyte.concentration_elyte(:, 1));
                
            battery.ThermodynamicParams.electrolyte.potential_mean_elyte_neg(1) = ...
                calculator.compute_electrolyte_mean_potential(battery, ...
                battery.ThermodynamicParams.electrolyte.potential_elyte(:, 1), "negative");

            battery.ThermodynamicParams.electrolyte.potential_mean_elyte_pos(1) = ...
                calculator.compute_electrolyte_mean_potential(battery, ...
                battery.ThermodynamicParams.electrolyte.potential_elyte(:, 1), "positive");
            % --------------------------------------------------------------------
            % 2. Main electrode loop
            % --------------------------------------------------------------------
            % -------- Electrode domain --------
            for electrode = working_electrode
                % ---- Prefix for electrodes ('neg', 'pos') ----------------------
                prefix = electrode{1}(1:3);
                Conc = battery.ConcentrationParams.electrode.(electrode);

                % -------- On aging resume: clear once and reuse last potential if needed
                if isLaterAging
                    Thermod = battery.ThermodynamicParams.electrode.(electrode);
                    Thermod.(['stoichiometry_' prefix]) = [];
                    Thermod.(['equilibrium_potential_' prefix]) = [];
                    Thermod.(['overpotential_' prefix]) = [];
                    battery.ThermodynamicParams.electrode.(electrode) = Thermod;
                end
                % ---- Stoichiometry -----------------------------------------------
                battery.ThermodynamicParams.electrode.(electrode).(['stoichiometry_' prefix])(1) = ...
                    calculator.compute_stoichiometry( ...
                    Conc.(['surface_concentration_' prefix])(1), ...
                    Conc.(['concentration_maximum_' prefix]));
                
                % ---- Equilibrium potential ---------------------------------------
                stoichiometry_init = battery.ThermodynamicParams.electrode.(electrode).(['stoichiometry_' prefix])(1);
                battery.ThermodynamicParams.electrode.(electrode).(['equilibrium_potential_' prefix])(1) = ...
                    calculator.compute_electrode_equilibrium_potential( ...
                        battery.ThermodynamicParams.electrode.(electrode).(['equilibrium_potential_equation_' prefix]), ...
                        stoichiometry_init);

                % ---- Aging stoichiometry ----------------------------------------
                battery.ThermodynamicParams.electrode.(electrode).(['aging_stoichiometry_' prefix])(1) = ...
                    calculator.compute_aging_stoichiometry( ...
                    Conc.(['average_concentration_' prefix])(1), ...
                    Conc.(['concentration_maximum_' prefix]));

                % ---- SEI resistance (first vs last for aging) ------------------
                if isLaterAging
                    sei_film_resistance = battery.ElectricalParams.electrode.(electrode).(['sei_film_resistance_' prefix])(end);
                else
                    sei_film_resistance = battery.ElectricalParams.electrode.(electrode).(['sei_film_resistance_' prefix])(1);
                end

                % ---- Overpotential --------------------------------------------
                battery.ThermodynamicParams.electrode.(electrode).(['overpotential_' prefix])(1) = ...
                    calculator.compute_electrode_overpotential(battery, ...
                    battery.KineticParams.electrode.(electrode).(['molar_ionic_flux_' prefix])(1), ...
                    battery.KineticParams.electrode.(electrode).(['exchange_current_density_' prefix])(1), ...
                    battery.ThermodynamicParams.electrolyte.(['potential_mean_elyte_' prefix])(1), ...
                    sei_film_resistance, ...
                    electrode);
                % --- Initial potential -----------------------------------------
                if isLaterAging
                    potential_last = battery.ThermodynamicParams.electrode.(electrode).(['potential_' prefix])(end);
                    battery.ThermodynamicParams.electrode.(electrode).(['potential_' prefix]) = [];
                    battery.ThermodynamicParams.electrode.(electrode).(['potential_' prefix])(1) = potential_last;
                else
                    battery.ThermodynamicParams.electrode.(electrode).(['potential_' prefix])(1) = ...
                        battery.ThermodynamicParams.electrode.(electrode).(['equilibrium_potential_' prefix])(1);
                end
            end
        end
    end
end