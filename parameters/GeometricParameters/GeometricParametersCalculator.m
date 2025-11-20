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
    end
end