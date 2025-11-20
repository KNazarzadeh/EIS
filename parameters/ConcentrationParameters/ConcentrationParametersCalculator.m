classdef ConcentrationParametersCalculator
    % Compute dependent battery parameters
    %
    % Description:
    % -----------
    % This class computes dependent parameters for batteries
    %
    % Properties:
    % -----------
    % params: Structure containing battery parameters
            
    methods
        function obj = ConcentrationParametersCalculator()

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     ELECTRODE FUNCTIONS     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% ------------------------ Electrode Concentration At Spcific Stoichiometry Level -------------------------- %%
        function concentration_at_sto = compute_electrode_concentration_at_stoichiometry_level(~, battery, electrode, stoichiometry_level)

            % Prefix for electrodes ('neg', 'pos')
            prefix = electrode{1}(1:3);

            stoichiometry_level = char(stoichiometry_level);

            Conc = battery.ConcentrationParams.electrode.(electrode);

            stoichiometry = battery.ThermodynamicParams.electrode.(electrode).(['stoichiometry_at_', stoichiometry_level, '_soc_', prefix])(end);

            concentraion_max = Conc.(['concentration_maximum_' prefix]);

            concentration_at_sto = stoichiometry * concentraion_max;
                    
        end

        %% ------------------------ Electrode Initial Concentration -------------------------- %%
        function electrod_concentration_init = compute_electrode_initial_concentration(~, battery, ...
                electrode)

            SOC_init = battery.SimulationParams.SOC_initial;
            % Prefix for electrodes ('neg', 'pos')
            prefix = electrode{1}(1:3);

            Conc = battery.ConcentrationParams.electrode.(electrode);
            Geom = battery.GeometricParams.electrode.(electrode).particles;

            concentration_s0 = Conc.(['concentration_stoichiometry_at_0_soc_' prefix])(1);

            concentration_s100 = Conc.(['concentration_stoichiometry_at_100_soc_' prefix])(1);

            numRadialNodes = Geom.(['number_radial_nodes_' prefix]);

            concentration_range = concentration_s100 - concentration_s0;

            electrod_concentration_init(1: numRadialNodes, 1) = concentration_range * SOC_init + concentration_s0;

        end

        %% ------------------------ Electrode Bulk Concentration -------------------------- %%
        function concentration = compute_electrode_concentration(~, ...
                Acs_hat, ...
                boundaryConcentrationDiscreteTimeMatrix, ...
                concentration_previous, ...
                molar_ionic_flux)

            concentration = -(Acs_hat \ ...
                (boundaryConcentrationDiscreteTimeMatrix * (molar_ionic_flux)' + concentration_previous));

        end

        %% ------------------------ Electrode Surface Concentration -------------------------- %%
        function concentration_surface = compute_electrode_surface_concentration(~, ...
                concentration, ...
                Acs_bar)

            concentration_surface = Acs_bar * concentration;

        end
        
        %% ------------------------ Electrode Average Concentration -------------------------- %%
        function average_concentration = compute_electrode_average_concentration(~, battery, ...
                concentration, ...
                electrode)

            % Prefix for electrodes ('neg', 'pos')
            prefix = electrode{1}(1:3);

            Geom = battery.GeometricParams.electrode.(electrode).particles;

            particle_radius = Geom.(['particle_radius_' prefix]);
            numRadialNodes = Geom.(['number_radial_nodes_' prefix]);
            radialGridSpacing = Geom.(['radialGridSpacing_' prefix]);

            radialGridSpacingVec  = (0:numRadialNodes-1).' * radialGridSpacing;       % column vector of radial coordinates
            
            % trapezoidal integration of c(r)*r^2 over [0, R]
            integrand = concentration .* (radialGridSpacingVec.^2);
            w = trapz(radialGridSpacingVec, integrand);   % ≈ ∫ c(r) r^2 dr

            average_concentration = (3 / particle_radius^3) * w;

        end

        %##################################################################################################%
        %################################   ---  ELECTROLYTE FUNCTIONS  ---    ############################%
        %##################################################################################################%
        %% ------------------------ Electrolyte Concentration --------------------------- %%
        function concentration_elyte = compute_electrolyte_concentration(~, ...
                surfaceDiffusionElyteDiscreteTime, ...
                boundaryFluxElectrolyteDiscreteTime, ...
                molar_ionic_flux_neg, ...
                molar_ionic_flux_pos, ...
                Concentration_elyte_previous)            
                
            % Update concentration
            concentration_elyte = -(surfaceDiffusionElyteDiscreteTime \ ...
                (boundaryFluxElectrolyteDiscreteTime * [molar_ionic_flux_neg; molar_ionic_flux_pos] + ...
                Concentration_elyte_previous));
            
        end

        %% ------------------------ Mean Electrolyte Concentration --------------------------- %%
        function concentration_mean_elyte = compute_electrolyte_mean_concentration(~, battery, electrode, concentration_elyte)
            
            numSpatialNodes_neg = battery.GeometricParams.electrode.negative.number_spatial_nodes_neg;

            % Determine spatial node indices based on electrode type
            if strcmpi(electrode, "negative")
                idx_range = 1:numSpatialNodes_neg;

            elseif strcmpi(electrode, "positive")
                numSpatialNodes_sep = battery.GeometricParams.separator.number_spatial_nodes_sep;
                idx_range = (numSpatialNodes_neg + numSpatialNodes_sep + 1):size(concentration_elyte, 1);
            end

            % Compute mean concentration
            concentration_mean_elyte = mean(concentration_elyte(idx_range), 'all')';
        end

    end
end