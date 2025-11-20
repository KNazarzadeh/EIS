classdef ThermodynamicParametersCalculator
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

        function obj = ThermodynamicParametersCalculator()

        end

        %% ----------------------------- Compute Stoichiometries ---------------------------- %%
        function stoichiometry = compute_stoichiometry(~, surface_concentration, max_concentration)
                                               
            stoichiometry = surface_concentration ./ max_concentration;
        end

        %% -------------------- Electrode Equilibrium Potential Curve -------------------- %%
        function equilibrium_potential = compute_electrode_equilibrium_potential(~, ...
            equilibrium_potential_equation, ...
            stoichiometry)

            equilibrium_potential = equilibrium_potential_equation(stoichiometry);
            
        end

        %% -------------------- Solid-phase Potential ------------------ %%
        function potential = compute_electrode_potential(~, overpotential, equilibrium_potential)     

            potential = overpotential + sum(equilibrium_potential, 2);

        end  

        %% -------------------- Solid-phase OverPotential ------------------ %%
        function overpotential = compute_electrode_overpotential(~, battery, ...
                molar_ionic_flux, ...
                exchange_current_density, ...
                electrolyte_potential_mean, ...
                film_resistance, ...
                temperature)     

            if nargin < 8
                temperature =  battery.ThermalParams.cell.temperature_reference_cell;
            end

            constants = ConstantParameters();
            
            WeightedFluxExchange = molar_ionic_flux ./ exchange_current_density;

            overpotential = asinh(WeightedFluxExchange .* constants.F) .* ...
                            (constants.R_gas * temperature / constants.F) + ...
                            electrolyte_potential_mean + ...
                            constants.F .* film_resistance .* molar_ionic_flux;
        end

        %% ------------------------ Electrolyte Potential --------------------------- %%
        function potential_elyte = compute_electrolyte_potential(~, ...
                                                                 boundaryPotentialElyteDiscreteTime, ...
                                                                 surfacePotentialElyteMatrix, ...
                                                                 diffusionPotentialElyteMatrix, ...
                                                                 molar_ionic_flux_neg, ...
                                                                 molar_ionic_flux_pos, ...
                                                                 concentration_elyte_previous)

                potential_elyte = -(surfacePotentialElyteMatrix \ (boundaryPotentialElyteDiscreteTime * ...
                        [molar_ionic_flux_neg; molar_ionic_flux_pos] + ...
                        diffusionPotentialElyteMatrix * log(concentration_elyte_previous)));
        end

        %% ------------------------ Mean Electrolyte Potential --------------------------- %%
        function potential_mean_elyte = compute_electrolyte_mean_potential(~, battery, ...
                potential_elyte, ...
                electrode)

            Geom = battery.GeometricParams;
            % Determine spatial node indices based on electrode type
            numSpatialNodes_neg = Geom.electrode.negative.number_spatial_nodes_neg;
            numSpatialNodes_sep = Geom.separator.number_spatial_nodes_sep;

            if strcmpi(electrode, "negative")
                idx_range = 1:numSpatialNodes_neg;
            elseif strcmpi(electrode, "positive")
                idx_range = (numSpatialNodes_neg + numSpatialNodes_sep +1):size(potential_elyte, 1);
            end

            % Compute mean concentration
            potential_mean_elyte = mean(potential_elyte(idx_range), 'all')';

        end

        %% ----------------------------- Compute Aging Stoichiometries ---------------------------- %%
        function aging_stoichiometry = compute_aging_stoichiometry(~, average_concentration, ...
                max_concentration)                                  
            
            aging_stoichiometry = average_concentration ./ max_concentration;
        end

        %% ------------------------ Solid-phase Potential Time Derivation --------------------------- %%
        function potential_timederivation = compute_electrode_potential_derivation(~, potential_current, potential_prev, timeStep)

            potential_timederivation = (potential_current - potential_prev) / timeStep;

        end

        %% -------------- Solid-phase Overpotential Derivation by Molar Flux--------------------- %%
        function overpotential_molarderivation = compute_electrode_overpotential_molarderivation(~, ...
                battery, molar_ionic_flux, exchange_current_density, film_resistance)

            temperature = battery.ThermalParams.cell.temperature_reference_cell;

            constants = ConstantParameters();

            overpotential_molarderivation = (constants.R_gas * temperature / constants.F) * ...
                          (constants.F ./ exchange_current_density) / ...
                          sqrt(1 + (constants.F .* molar_ionic_flux ./ exchange_current_density).^2) - ...
                          constants.F * film_resistance;
        end

    end
end