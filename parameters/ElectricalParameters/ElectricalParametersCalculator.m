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
        function soc_init = compute_initial_SOC(~, battery, cycleNum)
            
            if cycleNum > 1
                sto0_neg = battery.ThermodynamicParams.electrode.negative.stoichiometry_at_0_soc_neg(end);
                sto100_neg = battery.ThermodynamicParams.electrode.negative.stoichiometry_at_100_soc_neg(end);
                concentraion_max_neg = battery.ConcentrationParams.electrode.negative.concentration_maximum_neg;
                cs_avg_neg = battery.ConcentrationParams.electrode.negative.average_concentration_neg(1);
                soc_init = (cs_avg_neg/concentraion_max_neg - sto0_neg)/(sto100_neg - sto0_neg);
            else
                soc_init = battery.SimulationParams.('SOC_initial');
                if soc_init > 2
                    z = linspace(0,1,10000);
                    OCV =  battery.ThermodynamicParams.electrode.(electrode).U_pos((p.s100_pos-p.s0_pos)*z + p.s0_pos) - p.U_neg((p.s100_neg-p.s0_neg)*z + p.s0_neg);
                    soc_init = interp1(OCV,z,V_init,'linear','extrap');
                else
                    soc_init = max(min(soc_init, 1-1e-6), 1e-6);   % e.g., 1e-6 instead of 0
                end
            end
        end

        %% --------------------  State of Charge (SOC) -------------------- %%
        function soc = compute_SOC(~, battery, ...
                applied_current, ...
                soc_previous, ...
                battery_capacity_cell, ...
                cycleNum, ...
                timeStep, ...
                time)

            if cycleNum > 1
                sto0_neg = battery.ThermodynamicParams.electrode.negative.stoichiometry_at_0_soc_neg(end);
                sto100_neg = battery.ThermodynamicParams.electrode.negative.stoichiometry_at_100_soc_neg(end);
                concentraion_max_neg = battery.ConcentrationParams.electrode.negative.concentration_maximum_neg;
                cs_avg_neg = battery.ConcentrationParams.electrode.negative.average_concentration_neg(time+1);
                soc = (cs_avg_neg/concentraion_max_neg - sto0_neg)/(sto100_neg - sto0_neg);
            else
                soc = soc_previous + timeStep * applied_current / battery_capacity_cell;
            end
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
        function battery_capacity = compute_battery_capacity(~, battery, cycleNum)

            constants = ConstantParameters();
            if cycleNum > 1
                battery_capacity = (battery.ThermodynamicParams.electrode.negative.stoichiometry_at_100_soc_neg - ...
                    battery.ThermodynamicParams.electrode.negative.stoichiometry_at_0_soc_neg) * ...
                    battery.GeometricParams.electrode.negative.active_material_volume_fraction_neg * ...
                    battery.GeometricParams.electrode.negative.thickness_neg * ...
                    battery.GeometricParams.cell.surface_area_cell * ...
                    constants.F * ...
                    battery.ConcentrationParams.electrode.negative.concentration_maximum_neg;
            else
                %Reversible capacity of the battery [C]
                battery_capacity = (battery.ThermodynamicParams.electrode.positive.stoichiometry_at_0_soc_pos - ...
                    battery.ThermodynamicParams.electrode.positive.stoichiometry_at_100_soc_pos) * ...
                    battery.GeometricParams.electrode.positive.active_material_volume_fraction_pos * ...
                    battery.GeometricParams.electrode.positive.thickness_pos * ...
                    battery.GeometricParams.cell.surface_area_cell * ...
                    constants.F * ...
                    battery.ConcentrationParams.electrode.positive.concentration_maximum_pos;
            end
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

        %% ------------------------ Cell State of Health (SOH) --------------------------- %%

        function SOH = compute_SOH(~, battery, aging_capacity_battery)
            
            SOH =  aging_capacity_battery / battery.ElectricalParams.cell.('battery_capacity_cell');
            
        end

        %% -------------- AGING CALCULATIONS -------------------------------------------------------
        % Capacity Loss by LAM (Loss Active Material)
        function capacityLossLAM = calculate_capacityLoss_LAM(~, battery, ...
                capacityLossLAM_previous, ...
                volumeFractionLAM_derivation, ...
                stoichemiotry, ...
                electrode, timeStep)

            % calculate_capacityLoss_LAMNeg
            % -------------------------------------------------------------------------
            % Calculates the instantaneous capacity loss due to Loss of Active Material
            % (LAM) in a lithium-ion battery Negative electrode.
            %
            % Equation:
            %   dQ_LAM/dt = (dε_s/dt) * stoichiometry_neg * V_porous_neg * c_s,max
            %
            % Inputs:
            %   volumeFractionLAM_derivation - vector [1/s]
            %       The rate of change of active material volume fraction (ε_s) at each 
            %       timestep. Typically computed as: dε_s/dt = K_LAM * |I|
            %
            %   stoichiometry - scalar [-]
            %       The stoichiometric coefficient (y) representing the fraction of 
            %       active sites contributing to lithium storage.
            %
            %   porousVolume - scalar [m^3]
            %       Volume of the porous electrode.
            %
            %   concentration_max - scalar [mol/m^3]
            %       Maximum lithium concentration in the active material
            %
            % Output:
            %   capacityLossLAM - vector [mol/s]
            %       Rate of capacity loss caused by LAM at each timestep.
            %
            % -------------------------------------------------------------------------
               
            prefix = electrode{1}(1:3);

            concenteration_max = battery.ConcentrationParams.electrode.(electrode).(['concentration_maximum_' prefix]);

            porousVolume = battery.GeometricParams.electrode.(electrode).(['porousVolume_' prefix]);

            % Compute LAM capacity loss at each timestep
            % --- Positive loss rate in mol/s: dQ_loss/dt = -(dε/dt)*Δθ*V*c_s,max ---
            capacityLossLAM = capacityLossLAM_previous + ...
                (volumeFractionLAM_derivation * stoichemiotry * porousVolume * concenteration_max) * timeStep;
            
        end

        %% Capacity Loss by SEI
        function capacityLossSEI = calculate_capacityLoss_SEI(~, battery, ...
                Lithium_loss_SEI_formation)
  
            porousVolume = battery.GeometricParams.electrode.negative.porousVolume_neg;

            capacityLossSEI = Lithium_loss_SEI_formation * porousVolume;
        end


        %% Film Resistance
        function filmResistance = compute_filmResistance(~, battery, ...
                SEI_inner_thickness_total, ...
                electrode ...
                )
   
            prefix = electrode{1}(1:3);

            outer_layer_thickness = battery.GeometricParams.electrode.(electrode).(['sei_outer_thickness_initial_' prefix]);
            
            SEI_total_thickness = outer_layer_thickness + SEI_inner_thickness_total;
            
            filmResistance = SEI_total_thickness / ...
                    battery.TransportParams.electrode.(electrode).(['sei_conductivity_' prefix]);
        end

        %% Total Capacity Loss = SEI + LAM loss
        function capacity_loss = compute_total_capacity_loss(~, ...
            capacityLossLAM, ...
            capacityLossSEI)

            capacity_loss = abs(capacityLossLAM) + capacityLossSEI;
        end

        %% Total Capacity Loss
        function total_capacityLoss = compute_total_capacityLoss(~, ...
                capacityLossLAM_neg, ...
                capacityLossLAM_pos, ...
                capacityLossSEI_neg)
            total_capacityLoss = abs(capacityLossLAM_neg + capacityLossLAM_pos) + capacityLossSEI_neg;
        end
    end
end