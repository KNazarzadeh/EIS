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
                applied_current, ...
                working_electrode, ...
                cycleNum)

            % INITIALGEOMETRICPARAMETERS  Populate GeometricParams efficiently
            %
            %   - All heavy struct look-ups are cached locally.
            %   - Particle-field names are built once (vectorised where possible).
            %   - Aging adjustments are performed in a single pass.
            %   - Model-specific geometry is added after the electrode loop.

            % --------------------------------------------------------------------
            % 1. Initialise calculator & aging flags
            % --------------------------------------------------------------------
            isFirstAging = (cycleNum == 1);
            isLaterAging = (cycleNum > 1);
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

                % ---- Aging parameters -----------------------------------------------------------
                if isFirstAging && (strcmpi(electrode, "negative"))
                    % ---- SEI initialization (negative only) -------------------------------------
                        Geom.SEI_inner_thickness_total_neg(1) = Geom.sei_inner_thickness_initial_neg;
                        Geom.SEI_inner_thickness_increased_neg(1) = 0;

                elseif isLaterAging
                    % ---- Active material volume fraction
                    active_material_volume_fraction = Geom.(['active_material_volume_fraction_' prefix])(end);
                    Geom.(['active_material_volume_fraction_' prefix]) = [];
                    Geom.(['active_material_volume_fraction_' prefix])(1) = active_material_volume_fraction;

                    % ---- Clear old fields -----------------------------
                    % ---- Specific interfacial surface area --------------------------------------
                    Geom.particles.(['specific_interfacial_surface_area_' prefix]) = [];
                    % ---- LAM derivation active material volume fraction -------------------------
                    Geom.(['volumeFractionLAM_derivation_' prefix]) = [];

                    % ---- SEI initialization (negative only) -------------------------------------
                    if strcmpi(electrode,'negative')
                        % inner thickness
                        sei_inner_thickness_initial = Geom.SEI_inner_thickness_total_neg(end);
                        Geom.sei_inner_thickness_initial_neg = sei_inner_thickness_initial;
                        Geom.SEI_inner_thickness_total_neg = [];
                        Geom.SEI_inner_thickness_total_neg(1) = sei_inner_thickness_initial;
                        Geom.SEI_inner_thickness_increased_neg = [];
                        Geom.SEI_inner_thickness_increased_neg(1) = 0;
                    end
                end

                % ---- LAM derivation (scalar, once per electrode) --------------------------------
                Geom.(['volumeFractionLAM_derivation_' prefix])(1) = ...
                    calculator.compute_electrode_volumeFractionLAM_derivation(battery, ...
                    applied_current(1), ...
                    electrode);

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

                % ---- Porous volume per electrode ------------------------------------------------
                % specific_interfacial_surface_area = Geom.particles.(['surface_specific_interfacial_surface_area_' prefix]);            
                Geom.(['porousVolume_' prefix]) = ...
                    calculator.compute_electrode_porous_volume(battery, ...
                    specific_interfacial_surface_area, ...
                    thickness);

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