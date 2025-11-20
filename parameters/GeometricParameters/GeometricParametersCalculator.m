classdef GeometricParametersCalculator
    % Compute dependent battery parameters
    %
    % Description:
    % -----------
    % This class computes dependent parameters for batteries 
    methods
        function obj = GeometricParametersCalculator()
        end

        %% ------------------------ Electrode Particle Radial GridSpacing ----------------------- %%
        function radial_gridSpacing = compute_particle_radial_gridSpacing(~, ...
                numRadialNodes, ...
                particleRadius)
            
            radial_gridSpacing = particleRadius ./ (numRadialNodes - 1);
        end

        %% ----------------------------- Electrode Particle Surface ----------------------------- %%
        function surface_particle = compute_surface_particle(~, battery, electrode)

            % Prefix for electrodes ('neg', 'pos')
            prefix = electrode{1}(1:3);

            particleRadiusVec = battery.GeometricParams.electrode.(electrode).particles.(['particle_radius_' prefix]);

            particleFractionVec = 1;

            surface_particle = sum(particleFractionVec .* particleRadiusVec)/sum(particleFractionVec);
        end

        %% -------------------------- Specific Interfacial Surface Area ------------------------- %%
        function specific_interfacial_surface_area = compute_specific_interfacial_surface_area(~, ...
                particleRadius, ...
                active_material_volume_fraction)

            specific_interfacial_surface_area = 3 * active_material_volume_fraction ./ particleRadius;
        end

        %% -------------------------- Surface Specific Interfacial Surface Area ------------------------- %%
        function surface_specific_interfacial_surface_area = ...
                compute_surface_specific_interfacial_surface_area(~, ...
                    surface_particle_radius, ...
                    active_material_volume_fraction)
            
            surface_specific_interfacial_surface_area = 3 * active_material_volume_fraction ./ ...
                surface_particle_radius;
        end

        %% -------------------------------- Spatial GridSpacing --------------------------------- %%
        function spatial_gridSpacing = compute_spatial_gridSpacing(~, thickness, numSpatialNodes)
            
            spatial_gridSpacing = thickness / numSpatialNodes;

        end
        %% ------------------------ Electrolyte Volume Fraction ------------------------%%
        function electrolyteVolFraction = compute_electrolyteVolFraction(~, ...
                active_material_volume_fraction, ...
                filler_binder_volume_fraction)

            electrolyteVolFraction = 1 - ...
                (active_material_volume_fraction/ (1 - filler_binder_volume_fraction));
        end

        %% -------------- Electrode Volume -------------------------------------------------------
        function compute_electrode_volume = compute_electrode_volume(~, cell_surface_area, thickness)
           
            compute_electrode_volume = cell_surface_area * thickness;
            
        end     

        %% -------------- Electrode Surface Area -------------------------------------------------------
        function compute_electrode_area = compute_electrode_surface_area(~, battery, electrode)
            
            % Prefix for electrodes ('neg', 'pos', 'sep')
            prefix = electrode{1}(1:3);
            
            electrode_shape = battery.GeometricParams.electrode.(electrode).(['electrode_shape_'  prefix]);

            % Calculate area based on shape
            if strcmpi(electrode_shape, "rectangular") || strcmpi(electrode_shape, "regtangle")

                compute_electrode_area = battery.GeometricParams.electrode.(electrode).(['width' prefix]) * ...
                    battery.GeometricParams.electrode.(electrode).(['length_' prefix]);

            elseif strcmpi(electrode_shape, "circle") || strcmpi(electrode_shape, "disc")

                compute_electrode_area = pi * ...
                    (battery.GeometricParams.electrode.(electrode).(['diameter_' prefix])/2) ^ 2;
            end
        end    

        %% ================================ AGING CALCULATIONS =====================================
        % ------------------------------------ Electrode Porous Volume -----------------------------
        function electrode_porous_volume = compute_electrode_porous_volume(~, battery, ...
                specific_interfacial_surface_area, ...
                thickness)
           
            cell_surface_area = battery.GeometricParams.cell.('surface_area_cell');
            
            % Compute derivative of active material volume fraction at each timestep
            electrode_porous_volume = cell_surface_area * thickness * specific_interfacial_surface_area;
        end

        %% Electrode active material volume fraction due to Loss of Active Material (LAM)
        function volumeFraction = compute_electrode_volumeFraction(~, ...
            active_material_volume_fraction_previous, ...
            volumeFractionLAM_derivation, ...
            timeStep)

            % Compute derivative of active material volume fraction at each timestep
            volumeFraction = active_material_volume_fraction_previous + ...
                volumeFractionLAM_derivation * timeStep;
        end

        %% derivative of the electrode active material volume fraction due to Loss of Active Material (LAM)
        function volumeFractionLAM_derivation = compute_electrode_volumeFractionLAM_derivation(~, ...
            battery, ...
            applied_current, ...
            electrode ...
            )
            % compute_electrode_volumeFractionLAM_derivation
            % -------------------------------------------------------------------------
            % Calculates the time derivative of the electrode active material volume 
            % fraction due to Loss of Active Material (LAM).
            %
            % Equation:
            %   d(ε_s)/dt = -LAM_rate_loss * |I|
            %
            % Inputs:
            %   applied_current - vector [A] 
            %       Applied current profile of the cell (positive for charge,
            %       negative for discharge). Each element corresponds to one timestep.
            %
            %   LAM_rate_loss   - scalar [1/A·s] 
            %       Proportionality constant representing the LAM degradation rate.
            %
            % Output:
            %   volumeFractionLAM_derivation - vector [1/s]
            %       The rate of change of active material volume fraction (ε_s) at each
            %       timestep, same size as applied_current.
            %
            % Example:
            %   I = [1, -2, 3];               % applied current [A]
            %   K_LAM = 1e-6;                 % LAM rate constant
            %   d_eps = compute_volumeFractionLAM_derivation(I, K_LAM);
            %
            % -------------------------------------------------------------------------
            
            % Prefix for electrodes ('neg', 'pos')
            prefix = electrode{1}(1:3);

            kLAM = battery.KineticParams.electrode.(electrode).(['lam_rate_loss_' prefix]);
            % Compute derivative of active material volume fraction at each timestep
            volumeFractionLAM_derivation = -kLAM .* abs(applied_current);
        end

        %% SEI Increased Thickness
        function SEI_inner_thickness_increased = compute_SEI_inner_thickness_increased(~, battery, ...
                SEI_tunneling_reaction_rate, ...
                SEI_thickness_ratio, ...
                timeStep)
            
            Geom = battery.GeometricParams.electrode.negative;

            sei_molar_mass = Geom.sei_molar_mass_neg;
            sei_density = Geom.sei_density_neg;

            SEI_inner_thickness_increased = timeStep * ...
                (SEI_thickness_ratio * SEI_tunneling_reaction_rate * sei_molar_mass) / ...
                (2 * sei_density);
        end

        %% SEI Total Thickness
        function SEI_inner_thickness_total = compute_SEI_inner_thickness_total(~, ...
                SEI_inner_thickness_total_previous, ...
                SEI_inner_thickness_increased)


            SEI_inner_thickness_total = SEI_inner_thickness_total_previous + ...
                SEI_inner_thickness_increased;
        end

        %% Percentage of inner SEI thickness on the total SEI thickness
        function SEI_inner_thickness_ratio = compute_SEI_inner_thickness_ratio(~, ...
                battery, ...
                SEI_thickness_total)

            sei_thickness_initial_neg = battery.GeometricParams.electrode.negative.sei_thickness_initial_neg;

            SEI_inner_thickness_ratio = SEI_thickness_total ./ sei_thickness_initial_neg;
        end
    end
end