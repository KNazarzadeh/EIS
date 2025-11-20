classdef GeometricParameters
    %Initial and organize geometric parameters for battery types
    %
    % Description:
    % -----------
    % This class extracts geometric parameters from a specified battery type
    % and organizes them into domains (cell, electrode, electrolyte, separator, sei)
    % based on the electrochemical model type
    % (e.g., SPM, ESPM, DFN).
    %
   methods(Static)
        function initialGeometricParameters( ...
                battery, ...
                working_electrode)

            % INITIALGEOMETRICPARAMETERS  Populate GeometricParams efficiently
            %
            %   - All heavy struct look-ups are cached locally.
            %   - Particle-field names are built once (vectorised where possible).
            %   - Model-specific geometry is added after the electrode loop.
            % --------------------------------------------------------------------
            % 1. Initialise calculator
            % --------------------------------------------------------------------
            calculator = GeometricParametersCalculator();
            % --------------------------------------------------------------------
            % 2. Main electrode loop
            % --------------------------------------------------------------------
            % -------- Electrode domain --------
            for electrode = working_electrode
                % ---- Prefix for electrodes ('neg', 'pos') ----------------------
                prefix = electrode{1}(1:3);

                Geom = battery.GeometricParams.electrode.(electrode); % <-- local cache
                % ---- Number of Radial Nodes ---------------------------
                numRadialNodes = Geom.particles.(['number_radial_nodes_' prefix]);

                numSpatialNodes = Geom.(['number_spatial_nodes_' prefix]);

                % ---- Electrode Thickness ------------------------------
                thickness = Geom.(['thickness_' prefix]);

                % ---- Extract all particle Radius in one vector --------
                particleRadius = Geom.particles.(['particle_radius_' prefix]);

                % ---- Specific interfacial surface area --------------------
                specific_interfacial_surface_area = calculator.compute_specific_interfacial_surface_area( ...
                        particleRadius, ...
                        Geom.(['active_material_volume_fraction_' prefix])(1));

                Geom.particles.(['specific_interfacial_surface_area_' prefix])(1) = ...
                    specific_interfacial_surface_area;

                % ---- Radial Grid Spacing ----------------------------------
                Geom.particles.(['radialGridSpacing_' prefix]) = ...
                     calculator.compute_particle_radial_gridSpacing(numRadialNodes, particleRadius);

                % ---- Spatial Grid Spacing -----------------------
                Geom.(['spatialGridSpacing_' prefix]) = ...
                    calculator.compute_spatial_gridSpacing( ...
                    thickness, ...
                    numSpatialNodes);

                % ---- Electrolyte volume fraction (fallback computation) -------------------------
                elyteVfld = ['electrolyte_volume_fraction_' prefix];
                if ~isfield(Geom, (elyteVfld)) || isnan(Geom.(elyteVfld))
                    Geom.(['electrolyte_volume_fraction_' prefix]) = ...
                        calculator.compute_electrolyteVolFraction( ...
                        Geom.(['active_material_volume_fraction_' prefix]), ...
                        Geom.(['filler_binder_volume_fraction_' prefix]));
                end

                % --------------------------------------------------------
                % 3. Write back the fully-populated local struct
                % --------------------------------------------------------
                battery.GeometricParams.electrode.(electrode) = Geom;
            end
            % --------------------------------------------------------------------
            % 4. Model-specific geometry (executed once, after electrodes)
            % --------------------------------------------------------------------
            % ESPM Model: Extended Single Particle Model (ESPM)
            % --- Separator: Spatial Grid Spacing -----------------------------------------
            thickness = battery.GeometricParams.separator.thickness_sep; 
            numSpatialNodes_sep = battery.GeometricParams.separator.number_spatial_nodes_sep;
            battery.GeometricParams.separator.spatialGridSpacing_sep = ...
            calculator.compute_spatial_gridSpacing(thickness, numSpatialNodes_sep);
        end
    end
end