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

        %%
        function [sto_0_neg, sto_100_neg, sto_0_pos, sto_100_pos] = update_stoichiometry_levels(~, battery, total_capacityLoss_last)

            % --- parameters ---
            constants = ConstantParameters();
            Vmax = battery.ElectricalParams.cell.upper_voltage_cell;
            Vmin = battery.ElectricalParams.cell.lower_voltage_cell;
        
            Un = battery.ThermodynamicParams.electrode.negative.equilibrium_potential_equation_neg;  % Un(theta_n)
            Up = battery.ThermodynamicParams.electrode.positive.equilibrium_potential_equation_pos;  % Up(theta_p)
        
            epss_n = battery.GeometricParams.electrode.negative.active_material_volume_fraction_neg;
            epss_p = battery.GeometricParams.electrode.positive.active_material_volume_fraction_pos;
            A      = battery.GeometricParams.cell.surface_area_cell;
        
            csmax_n = battery.ConcentrationParams.electrode.negative.concentration_maximum_neg;
            csmax_p = battery.ConcentrationParams.electrode.positive.concentration_maximum_pos;
            del_n   = battery.GeometricParams.electrode.negative.thickness_neg;
            del_p   = battery.GeometricParams.electrode.positive.thickness_pos;
        
            Qneg = mean(epss_n) * del_n * csmax_n * A * constants.F;
            Qpos = mean(epss_p) * del_p * csmax_p * A * constants.F;
        
            % use the cycle-1 100% sto values to define QLi0 (same as your previous code)
            s100n0 = battery.ThermodynamicParams.electrode.negative.stoichiometry_at_100_soc_neg(1,:);
            s100p0 = battery.ThermodynamicParams.electrode.positive.stoichiometry_at_100_soc_pos(1,:);
            QLi0   = s100n0*Qneg + s100p0*Qpos;
            QLi    = QLi0 - total_capacityLoss_last;
        
            % seeds from last known sto (keeps continuity if available)
            tn0_seed   = battery.ThermodynamicParams.electrode.negative.stoichiometry_at_0_soc_neg(1,:);
            tn100_seed = battery.ThermodynamicParams.electrode.negative.stoichiometry_at_100_soc_neg(1,:);
        
            % tiny interior box (avoid exact 0/1 for OCV eval)
            th_lo = 1e-9; th_hi = 1-1e-9;
        
            % feasible theta_n interval from tp in (th_lo, th_hi)
            tn_min = max(th_lo, (QLi - (1-1e-9)*Qpos)/Qneg);
            tn_max = min(th_hi, (QLi - 1e-9*Qpos)/Qneg);
        
            % if infeasible, project QLi to nearest feasible band and recompute interval
            if ~(tn_min < tn_max)
                QLi = min(max(QLi, th_lo*(Qneg+Qpos)), (1-1e-9)*(Qneg+Qpos));
                tn_min = max(th_lo, (QLi - (1-1e-9)*Qpos)/Qneg);
                tn_max = min(th_hi, (QLi - 1e-9*Qpos)/Qneg);
                if ~(tn_min < tn_max)
                    % pathological: return mid values to keep sim running
                    tn_mid = 0.5*(th_lo+th_hi);
                    tp_mid = (QLi - tn_mid*Qneg)/Qpos;
                    sto_0_neg = tn_mid;  sto_100_neg = tn_mid;
                    sto_0_pos = min(max(tp_mid, th_lo), th_hi);
                    sto_100_pos = sto_0_pos;
                    return;
                end
            end
        
            % helper lambdas (anonymous; no extra functions)
            g_at = @(V) @(tn) Up( (QLi - tn.*Qneg)./Qpos ) - Un(tn) - V;   % scalar residual g(tn)=0
            clamp = @(x) min(max(x, th_lo), th_hi);
        
            % ---- solve for 100% (Vmax) ----
            g100 = g_at(Vmax);
            % bracket on feasible interval
            g_lo = g100(tn_min); g_hi = g100(tn_max);
            opts = optimset('TolX',1e-12,'TolFun',1e-12,'Display','off');
        
            if sign(g_lo)*sign(g_hi) <= 0
                % root is inside feasible interval
                try
                    tn100 = fzero(g100, [tn_min, tn_max], opts);
                catch
                    % fallback: start near previous seed if available
                    tn100 = fzero(g100, clamp(tn100_seed), opts);
                end
            else
                % no root → Vmax outside achievable OCV span → pick nearest boundary
                tn100 = (abs(g_lo) <= abs(g_hi)) * tn_min + (abs(g_lo) > abs(g_hi)) * tn_max;
            end
            tn100 = clamp(tn100);
            tp100 = clamp( (QLi - tn100*Qneg)/Qpos );
        
            % ---- solve for 0% (Vmin) ----
            g0 = g_at(Vmin);
            g_lo = g0(tn_min); g_hi = g0(tn_max);
        
            if sign(g_lo)*sign(g_hi) <= 0
                try
                    tn0 = fzero(g0, [tn_min, tn_max], opts);
                catch
                    tn0 = fzero(g0, clamp(tn0_seed), opts);
                end
            else
                tn0 = (abs(g_lo) <= abs(g_hi)) * tn_min + (abs(g_lo) > abs(g_hi)) * tn_max;
            end
            tn0 = clamp(tn0);
            tp0 = clamp( (QLi - tn0*Qneg)/Qpos );
        
            % outputs
            sto_0_neg   = tn0;    sto_100_neg = tn100;
            sto_0_pos   = tp0;    sto_100_pos = tp100;
        
            % quick sanity
            if sto_0_neg >= sto_100_neg
                % not fatal; can happen if both projected to boundaries
                % warning('update_stoichiometry_levels:ordering','tn0 >= tn100; check QLi/V window.');
            end
        end

    end
end