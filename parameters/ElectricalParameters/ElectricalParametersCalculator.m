classdef ElectricalParametersCalculator
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

        function obj = ElectricalParametersCalculator()

        end
        
        %% -------------------- Set Initial State of Charge (SOC) -------------------- %%
        function soc_init = compute_initial_SOC(~, battery)
            
            soc_init = battery.SimulationParams.('SOC_initial');
            if soc_init > 2
                z = linspace(0,1,10000);
                OCV =  battery.ThermodynamicParams.electrode.(electrode).U_pos((p.s100_pos-p.s0_pos)*z + p.s0_pos) - p.U_neg((p.s100_neg-p.s0_neg)*z + p.s0_neg);
                soc_init = interp1(OCV,z,V_init,'linear','extrap');
            else
                soc_init = max(min(soc_init, 1-1e-6), 1e-6);   % e.g., 1e-6 instead of 0
            end
        end

        %% --------------------  State of Charge (SOC) -------------------- %%
        function soc = compute_SOC(~, ...
                applied_current, ...
                soc_previous, ...
                battery_capacity_cell, ...
                timeStep)

                soc = soc_previous + timeStep * applied_current / battery_capacity_cell;
        end

        %% ------------------------ Cell OCV Range Limitation --------------------------- %%
        function [OCV_min, OCV_max] = compute_ocv_range(~, battery)

            equilibrium_potential_equation_neg =  battery.ThermodynamicParams.electrode.negative.equilibrium_potential_equation_neg;
            equilibrium_potential_equation_pos =  battery.ThermodynamicParams.electrode.positive.equilibrium_potential_equation_pos;

            stoichiometry_at_0_soc_neg = battery.ThermodynamicParams.electrode.negative.stoichiometry_at_0_soc_neg;
            stoichiometry_at_100_soc_neg = battery.ThermodynamicParams.electrode.negative.stoichiometry_at_100_soc_neg;
            stoichiometry_at_0_soc_pos = battery.ThermodynamicParams.electrode.positive.stoichiometry_at_0_soc_pos;
            stoichiometry_at_100_soc_pos = battery.ThermodynamicParams.electrode.positive.stoichiometry_at_100_soc_pos;

            OCV_max = equilibrium_potential_equation_pos(stoichiometry_at_100_soc_pos) - equilibrium_potential_equation_neg(stoichiometry_at_100_soc_neg); 
            OCV_min = equilibrium_potential_equation_pos(stoichiometry_at_0_soc_pos) - equilibrium_potential_equation_neg(stoichiometry_at_0_soc_neg); 

        end

        %% ------------------------ Cell OCV Range Limitation --------------------------- %%
        function battery_capacity = compute_battery_capacity(~, battery)

            constants = ConstantParameters();
            %Reversible capacity of the battery [C]
            battery_capacity = (battery.ThermodynamicParams.electrode.positive.stoichiometry_at_0_soc_pos - ...
                battery.ThermodynamicParams.electrode.positive.stoichiometry_at_100_soc_pos) * ...
                battery.GeometricParams.electrode.positive.active_material_volume_fraction_pos * ...
                battery.GeometricParams.electrode.positive.thickness_pos * ...
                battery.GeometricParams.cell.surface_area_cell * ...
                constants.F * ...
                battery.ConcentrationParams.electrode.positive.concentration_maximum_pos;
        end
 
        %% ------------------------ Cell Voltage --------------------------- %%
        function voltage = compute_cell_voltage(~, battery, ...
                applied_current, time)
                            
            potential_pos = battery.ThermodynamicParams.electrode.positive.('potential_pos')(time+1);
            potential_neg = battery.ThermodynamicParams.electrode.negative.('potential_neg')(time+1);
           
            Rseries = (battery.ElectricalParams.electrode.positive.('contact_resistance_pos') + ...
                      battery.ElectricalParams.electrode.negative.('contact_resistance_neg')) / ...
                      battery.GeometricParams.cell.('surface_area_cell');

            voltage = potential_pos - potential_neg - applied_current * Rseries;

        end
        
        %% ------------------------ Cell Capacity --------------------------- %%
        function cell_capacity = compute_Cell_Capacity(~, ...
                SOC_previous, ...
                battery_capacity_cell, ...
                total_capacityLoss)

            % cell_capacity  = abs(battery.ThermodynamicParams.cell.SOC - battery.ThermodynamicParams.cell.SOC(1)) * ...
            %     (battery.ElectricalParams.cell.('battery_capacity_cell')(end) / 3600); 
            
            cell_capacity  = abs(SOC_previous * battery_capacity_cell - total_capacityLoss);
            
            % Ah
            % cell_capacity_by_current = cumtrapz(battery.SimulationParams.timeVector, ...
            %     abs(applied_current))/3600;

            % surface_area_cm2 = battery.GeometricParams.cell.surface_area_cell * 10000;
            % cell_capacity_mAh_cm2 = 1000 * cell_capacity / surface_area_cm2;

        end

        %% ------------------------ Cell Capacity -------------------------- %%
        function cell_nominal_capacity = compute_cell_nominal_capacity(~, battery, positive_electrode)
            % capacity of the battery [C]


            % Prefix for electrodes ('pos')
            prefix = char(extractBefore(positive_electrode,4));

            constants = ConstantParameters();

            if ~isfield( battery.GeometricParams.electrode.(positive_electrode), (['surface_area_' prefix])) || ...
                    isempty( battery.GeometricParams.electrode.(positive_electrode).(['surface_area_' prefix]))

                 battery.GeometricParams.electrode.(positive_electrode).(['surface_area_' prefix]) = compute_electrode_surface_area(obj, battery, positive_electrode);

            end

            cell_nominal_capacity = ( battery.ThermodynamicParams.electrode.(positive_electrode).(['stoichiometry_at_0_soc_' prefix]) - ...
                             battery.ThermodynamicParams.electrode.(positive_electrode).(['stoichiometry_at_100_soc_' prefix])) * ...
                             battery.GeometricParams.electrode.(positive_electrode).(['active_material_volume_fraction_' prefix]) * ...
                             battery.GeometricParams.electrode.(positive_electrode).(['thickness_' prefix]) * ...
                             battery.GeometricParams.electrode.(positive_electrode).(['surface_area_' prefix]) * ...
                             battery.ConcentrationParams.electrode.(positive_electrode).(['concentration_maximum_' prefix]) * ...
                            constants.F;
        end

    end
end