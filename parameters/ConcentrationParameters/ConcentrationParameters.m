classdef ConcentrationParameters
    % Extract and organize concentration parameters for battery models.
    %
    % Description:
    % -----------
    % This class extracts concentration parameters from a specified battery type
    % and organizes them into domains (cell, electrode,
    % electrolyte, separator) based on the electrochemical model type (e.g., SPM, ESPM, DFN).
    %
    methods(Static)
        % Initial Values of Parameters
        function initialConcentrationParameters( ...
                battery, ...
                preComputedMatrix, ...
                working_electrode)

            % Initial and organize concentration parameters by domain
            %
            % Description:
            % ------------
            % Extracts concentration parameters, organizing them into substructures:
            % cell, electrode (negative, positive, including particle parameters),
            % electrolyte, separator, and seei based on the battery type and
            % the specified electrochemical model (SPM, ESPM, etc.). 
            %
            % Outputs:
            % --------
            % obj: Updated object with concentrationParams structure populated
            
            % --------------------------------------------------------------------
            % 1. Initialise calculator
            % --------------------------------------------------------------------
            calculator = ConcentrationParametersCalculator();

            % --------------------------------------------------------------------
            % 2. Main electrode loop
            % --------------------------------------------------------------------
            % -------- Electrode domain --------
            for electrode = working_electrode
                % ---- Prefix for electrodes ('neg', 'pos') ----------------------
                prefix = electrode{1}(1:3);
                Conc = battery.ConcentrationParams.electrode.(electrode);  % <-- cache

                % ---- Stoichiometry concentrations (0% and 100% SOC) -----------------                
                Conc.(['concentration_stoichiometry_at_0_soc_' prefix]) = ...
                    calculator.compute_electrode_concentration_at_stoichiometry_level(battery, electrode, "0");
                Conc.(['concentration_stoichiometry_at_100_soc_' prefix]) = ...
                    calculator.compute_electrode_concentration_at_stoichiometry_level(battery, electrode, "100");

                % --------------------------------------------------------
                % --- Write back the fully-populated local struct
                % --------------------------------------------------------
                battery.ConcentrationParams.electrode.(electrode) = Conc;
                Conc = battery.ConcentrationParams.electrode.(electrode);  % <-- cache
                % ---- Compute initial values --------------------------
                Conc.(['concentration_' prefix])(:, 1) = ...
                    calculator.compute_electrode_initial_concentration(battery, ...
                    electrode);
            
                Acs_bar = preComputedMatrix.(electrode).Acs_bar;
                Conc.(['surface_concentration_' prefix])(1) = ...
                    calculator.compute_electrode_surface_concentration(...
                    Conc.(['concentration_' prefix])(:, 1), ...
                    Acs_bar);

                Conc.(['average_concentration_' prefix])(1) = ...
                    calculator.compute_electrode_average_concentration(battery, ...
                    Conc.(['concentration_' prefix])(:, 1), ...
                    electrode);
                % --------------------------------------------------------
                % 3. Write back the fully-populated local struct
                % --------------------------------------------------------
                battery.ConcentrationParams.electrode.(electrode) = Conc;
            end

            % --------------------------------------------------------------------
            % 4. Model-specific geometry (executed once, after electrodes)
            % --------------------------------------------------------------------
            % ESPM Model: Extended Single Particle Model (ESPM) uses electrolyte concentration
            Geom = battery.GeometricParams;
            ConcElyte = battery.ConcentrationParams.electrolyte;
            
            numSpatialNodes_sep = Geom.separator.('number_spatial_nodes_sep');
            % ---- Full Cell --------------------------------------------------------------
            % ---- Electrolyte Mean Concentration -------------------------------------
            numSpatialNodes_neg = Geom.electrode.negative.('number_spatial_nodes_neg');
            numSpatialNodes_pos = Geom.electrode.positive.('number_spatial_nodes_pos');
            % ---- Electrolyte parameters ---------------------------------------------------
            numSpatialNodesTotal = numSpatialNodes_neg + numSpatialNodes_sep + numSpatialNodes_pos;
            % ---- Electrolyte Concentration --------------------------------------
            concentration_elyte_init = ConcElyte.('concentration_initial_elyte') .* ones(numSpatialNodesTotal, 1);
            ConcElyte.('concentration_elyte')(:, 1) = concentration_elyte_init;
            ConcElyte.('concentration_initial_elyte') = concentration_elyte_init;
            
            idx_range_neg = 1:numSpatialNodes_neg;
            ConcElyte.('concentration_mean_elyte_neg')(1) = mean(concentration_elyte_init(idx_range_neg));
            idx_range_pos = (numSpatialNodes_neg + numSpatialNodes_sep + 1);
            ConcElyte.('concentration_mean_elyte_pos')(1) = mean(concentration_elyte_init(idx_range_pos, end));
            % --------------------------------------------------------
            % 5. Write back the fully-populated local struct
            % --------------------------------------------------------
            battery.GeometricParams = Geom;
            battery.ConcentrationParams.electrolyte = ConcElyte;
        end
    end    
end