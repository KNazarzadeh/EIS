classdef KineticParametersCalculator
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
        function obj = KineticParametersCalculator()
            % Empty constructor; battery will be passed directly
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     ELECTRODE FUNCTIONS     %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% ------------------------ Electrode Initial Molar Ionic Flux -------------------------- %%
        function molar_ionic_flux_init = compute_initial_molarIonicFlux(~, battery, ...
                                                                       applied_current, ...
                                                                       electrode)

            % Prefix for electrodes ('neg', 'pos')
            prefix = electrode{1}(1:3);

            % Compute the constant denominator
            constants = ConstantParameters();

            Geom = battery.GeometricParams.electrode.(electrode);
            
            surface_area = battery.GeometricParams.cell.('surface_area_cell');
            thickness = Geom.(['thickness_' prefix]);
            
            specific_interfacial_surface_area = Geom.particles.(['specific_interfacial_surface_area_' prefix]);
            denominator = surface_area * thickness * specific_interfacial_surface_area;
            
            % Calculate molar ionic flux for each electrode
            if strcmpi(electrode, "negative")
                molar_ionic_flux_init = -applied_current(1) / (constants.F * denominator);
    
            elseif strcmpi(electrode, "positive")
                molar_ionic_flux_init = applied_current(1) / (constants.F * denominator);
            end
        end

        %% -------------------------- Molar Ionic Flux Best Best ------------------------- %%
        function molar_ionic_flux = compute_molar_ionic_flux(~, battery, ...
                                                           applied_current, ...
                                                           specific_interfacial_surface_area, ...
                                                           electrode)

            % Prefix for electrodes ('neg', 'pos')
            prefix = electrode{1}(1:3);

            % Compute the constant denominator
            constants = ConstantParameters();
            Geom = battery.GeometricParams;

            surface_area = Geom.cell.('surface_area_cell');
            thickness = Geom.electrode.(electrode).(['thickness_' prefix]);

            denominator = surface_area * thickness * specific_interfacial_surface_area;

            % Calculate molar ionic flux for each electrode
            if strcmpi(electrode, "negative")
                molar_ionic_flux = -applied_current / (constants.F * denominator);
    
            elseif strcmpi(electrode, "positive")
                molar_ionic_flux = applied_current / (constants.F * denominator);
            end
        end

        %% -------------------------- Exchange Current Density ------------------------- %%
        function exchange_current_density = compute_exchange_current_density(~, battery, ...
                electrode_surface_concentration, ...
                concentration_mean_elyte, ...
                electrode)

            % Prefix for electrodes ('neg', 'pos')
            prefix = electrode{1}(1:3);
                        
            electrode_concentration_max = battery.ConcentrationParams.electrode.(electrode).(['concentration_maximum_' prefix]);
            
            kinet = battery.KineticParams.electrode;
            reaction_rate = kinet.(electrode).(['reaction_rate_coefficient_' prefix]);

            exchange_factor = (concentration_mean_elyte .^ kinet.negative.('charge_transfer_coefficient_neg')) .* ...
                ((electrode_concentration_max - electrode_surface_concentration) .^ kinet.negative.('charge_transfer_coefficient_neg')) .* ...
                (electrode_surface_concentration .^ kinet.positive.('charge_transfer_coefficient_pos'));

            exchange_current_density = reaction_rate .* exchange_factor;

        end

        %% -------------------------- Reaction Rate Coefficient ------------------------- %%
        function reaction_rate = compute_reaction_rate_coefficient(~, battery, electrode, temperature)

            temperature_ref = battery.ThermalParams.cell.('temperature_reference_cell');

            if nargin < 4
                temperature = temperature_ref;
            end

            % Prefix for electrodes ('neg', 'pos')
            prefix = electrode{1}(1:3);

            Kinet = battery.KineticParams.electrode.(electrode);

            reaction_rate_energy = Kinet.(['reaction_rate_activation_energy_' prefix]);
            reaction_rate_reference = Kinet.(['reaction_rate_coefficient_reference_' prefix]);

            if isnan(reaction_rate_energy) || isempty(reaction_rate_energy)
                reaction_rate = reaction_rate_reference;
                return;
            end

            % Compute Arrhenius factor
            arrhenius_factor = arrhenius_temperature_adjustment(reaction_rate_energy, ...
                                                                temperature_ref, ...
                                                                temperature);

            % if (isa(reaction_rate_reference,'function_handle'))
            %     % Get diffusion reference equation
            %     reaction_rate = @(ce, T) derivative_activity_reference_elyte(ce) .* arrhenius_factor(ce, T);
            % else
            %     reaction_rate = derivative_activity_reference_elyte * arrhenius_factor;
            % end            

            reaction_rate =  reaction_rate_reference .* arrhenius_factor;
        end

    end
end