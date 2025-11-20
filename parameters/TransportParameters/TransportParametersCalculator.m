classdef TransportParametersCalculator
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
        function obj = TransportParametersCalculator()
            % Empty constructor; battery will be passed directly
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     ELECTRODE FUNCTIONS     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% -------------------- Electrode Diffusion Coeffecient -------------------- %%
        function electrode_diffusion_coefficient = compute_electrode_diffusion_coefficient(obj, ...
                ThermalParams, ...
                TransportParams, ...
                electrode, ...
                stoichiometry, ...
                temperature)
            
            temperature_ref = ThermalParams.cell.('temperature_reference_cell');

            if nargin < 6
                temperature = temperature_ref;
            end
            
            prefix = char(extractBefore(electrode,4));

            Trans = TransportParams.electrode.(electrode);

            diffusion_activation_energy = Trans.(['diffusion_activation_energy_' prefix]);

            if ~isnan(diffusion_activation_energy)
                arrhenius_factor = arrhenius_temperature_adjustment(diffusion_activation_energy, ...
                                                                    temperature_ref, ...
                                                                    temperature);
            else
                arrhenius_factor = 1;
            end

            diffusion_coefficient_reference = Trans.(['diffusion_coefficient_reference_' prefix]);

            if (isa(diffusion_coefficient_reference,'function_handle'))
                % Get diffusion reference equation
                if nargin(diffusion_coefficient_reference) == 1
                    electrode_diffusion_coefficient = diffusion_coefficient_reference(stoichiometry) * arrhenius_factor;
                elseif nargin(diffusion_coefficient_reference) == 2
                    electrode_diffusion_coefficient = diffusion_coefficient_reference(stoichiometry, temperature) * arrhenius_factor;
                end
            else
                electrode_diffusion_coefficient = diffusion_coefficient_reference * arrhenius_factor;
            end            
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     ELECTROLYTE FUNCTIONS     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% ---------------------- Electrolyte Transference Number of Ions ----------------------- %%
        function transference_ion_number_elyte = compute_transference_ion_number_electrolyte(obj, battery)
           
            Trans_elyte = battery.TransportParams.electrolyte;
            if ~isfield(Trans_elyte, ('transference_ion_number_elyte')) || ...
                isnan(Trans_elyte.('transference_ion_number_elyte'))

                if ~isnan(Trans_elyte.('transference_ion_number_elyte_equation'))
                    
                    transference_ion_number_elyte_equation = Trans_elyte.('transference_ion_number_elyte_equation');
    
                    transference_ion_number_elyte = transference_ion_number_elyte_equation(concentration_initial_elyte);
                end
    
           else
                transference_ion_number_elyte = Trans_elyte.('transference_ion_number_elyte');
            end
        end  

        %% --------------------  Electrolyte Bruggeman Correction Factor ------------------- %%
        function brug_correction = compute_brugCorrection(~, battery)

            numSpatialNodes_sep = battery.GeometricParams.separator.('number_spatial_nodes_sep');
            electrolyte_volume_fraction_sep = battery.GeometricParams.separator.('electrolyte_volume_fraction_sep');
            bruggeman_coefficient_sep = battery.TransportParams.separator.('bruggeman_coefficient_sep');

            numSpatialNodes_neg = battery.GeometricParams.electrode.negative.('number_spatial_nodes_neg');
            numSpatialNodes_pos = battery.GeometricParams.electrode.positive.('number_spatial_nodes_pos');

            electrolyte_fraction = [battery.GeometricParams.electrode.negative.('electrolyte_volume_fraction_neg') * ones(numSpatialNodes_neg, 1);
                                    electrolyte_volume_fraction_sep * ones(numSpatialNodes_sep, 1);
                                    battery.GeometricParams.electrode.positive.('electrolyte_volume_fraction_pos') * ones(numSpatialNodes_pos, 1)
                                    ];

            bruggeman_coefficient = [battery.TransportParams.electrode.negative.('bruggeman_coefficient_neg') * ones(numSpatialNodes_neg, 1);
                                    bruggeman_coefficient_sep * ones(numSpatialNodes_sep, 1);
                                    battery.TransportParams.electrode.positive.('bruggeman_coefficient_pos') * ones(numSpatialNodes_pos, 1)
                                    ];

            brug_correction = electrolyte_fraction .^ bruggeman_coefficient;

        end
        %% --------------------  Electrolyte Effective Diffusion Coefficient ------------------- %%
        function effective_diffusion_elyte = compute_effective_diffusion_electrolyte(obj, ...
                ThermalParams, ... 
                TransportParams, ...
                concentration_elyte, ...
                temperature)
            
            batteryCell = ThermalParams.cell;

            batteryElyte = TransportParams.electrolyte;

            temperature_ref = batteryCell.('temperature_reference_cell');

            brug_correction_elyte = batteryElyte.('brug_correction_elyte');

            % -------------------- Default temperature --------------------
            if nargin < 5 || isempty(temperature)
                temperature = temperature_ref;
            end
            
            % -------------------- Evaluate Arrhenius factor --------------------
            diffusion_energy = TransportParams.electrolyte.('diffusion_activation_energy_elyte');
            
            if (isa(diffusion_energy, 'function_handle') || isnumeric(diffusion_energy)) && ...
                    ~isnan(diffusion_energy)
                % Compute Arrhenius factor
                arrhenius_factor = arrhenius_temperature_adjustment(diffusion_energy, ...
                                                                    temperature_ref, ...
                                                                    temperature, ...
                                                                    concentration_elyte);  

            elseif isnan(diffusion_energy)
                arrhenius_factor = 1;
            end

            % -------------------- Evaluate diffusion_coefficient_reference(concentration_elyte, temperature) --------------------
            % Extract parameters
            diffusion_coeff_ref = TransportParams.electrolyte.('diffusion_coefficient_reference_elyte');
        
            % Initialize diffusion coefficient
            if isnumeric(diffusion_coeff_ref)
                % Scalar case
                diffusion_coefficient = diffusion_coeff_ref;
            elseif isa(diffusion_coeff_ref, 'function_handle')
                % Function handle case: check number of input arguments
                nargin_inputs = nargin(diffusion_coeff_ref);
                if nargin_inputs == 1
                    diffusion_coefficient = diffusion_coeff_ref(concentration_elyte);
                elseif nargin_inputs == 2
                    diffusion_coefficient = diffusion_coeff_ref(concentration_elyte, temperature);
                else
                    error('diffusion_coefficient_reference_elyte expects 1 or 2 inputs, got %d', n_inputs);
                end
            else
                error('diffusion_coefficient_reference_elyte must be a scalar or function handle');
            end

            
            % -------------------- Effective diffusion (element-wise) --------------------
            effective_diffusion_elyte = diffusion_coefficient .* arrhenius_factor .* brug_correction_elyte;

        end        

        %% -------------------- Electrolyte Effective Conduction Coefficient  ------------------- %%
        function effective_conductivity_elyte = compute_effective_conductivity_electrolyte(~, ...
                ThermalParams, ... 
                TransportParams, ...
                concentration_elyte, ...
                temperature)

            batteryCell = ThermalParams.cell;

            batteryElyte = TransportParams.electrolyte;

            temperature_ref = batteryCell.('temperature_reference_cell');

            brug_correction_elyte = batteryElyte.('brug_correction_elyte');

            % -------------------- Default temperature --------------------
            if nargin < 5 || isempty(temperature)
                temperature = temperature_ref;
            end
            
            % -------------------- Evaluate Arrhenius factor --------------------
            conductivity_energy = batteryElyte.('conductivity_activation_energy_elyte');
            
            if (isa(conductivity_energy, 'function_handle') || isnumeric(conductivity_energy)) && ...
                    ~isnan(conductivity_energy)
                % Compute Arrhenius factor
                arrhenius_factor = arrhenius_temperature_adjustment(conductivity_energy, ...
                                                                    temperature_ref, ...
                                                                    temperature, ...
                                                                    concentration_elyte);  

            elseif isnan(conductivity_energy)
                arrhenius_factor = 1;
            end

            % -------------------- Evaluate diffusion_coefficient_reference(concentration_elyte, temperature) --------------------
            % Extract parameters
            conductivity_coeff_ref = batteryElyte.('conductivity_coefficient_reference_elyte');
        
            % Initialize diffusion coefficient
            if isnumeric(conductivity_coeff_ref)
                % Scalar case
                conductivity_coefficient = conductivity_coeff_ref;
            elseif isa(conductivity_coeff_ref, 'function_handle')
                % Function handle case: check number of input arguments
                nargin_inputs = nargin(conductivity_coeff_ref);
                if nargin_inputs == 1
                    conductivity_coefficient = conductivity_coeff_ref(concentration_elyte);
                elseif nargin_inputs == 2
                    conductivity_coefficient = conductivity_coeff_ref(concentration_elyte, temperature);
                else
                    error('conductivity_coefficient_reference_elyte expects 1 or 2 inputs, got %d', n_inputs);
                end
            else
                error('conductivity_coefficient_reference_elyte must be a scalar or function handle');
            end

            
            % -------------------- Effective diffusion (element-wise) --------------------
            effective_conductivity_elyte = conductivity_coefficient .* arrhenius_factor .* brug_correction_elyte;

        end  
        
        %% -------------------- Electrolyte derivative_activity_coefficient_concentraion Coeffecient -------------------- %%
        function derivative_coefficient_elyte = compute_electrolyte_derivative_activity_coefficient(~, ...
                ThermalParams, ... 
                TransportParams, ...
                concentration_elyte, ...
                temperature ...
                )

            batteryCell = ThermalParams.cell;

            batteryElyte = TransportParams.electrolyte;

            temperature_ref = batteryCell.('temperature_reference_cell');

            % -------------------- Default temperature --------------------
            if nargin < 5 || isempty(temperature)
                temperature = temperature_ref;
            end
            
            % Extract parameters
            derivative_energy = batteryElyte.('derivative_activity_coefficient_concentraion_energy_elyte');
            derivative_coeff_ref = batteryElyte.('derivative_activity_coefficient_concentraion_reference_elyte');

            if ~isa(derivative_energy, 'function_handle') && ...
                (isnan(derivative_energy) || isempty(derivative_energy))
                arrhenius_factor = 1;
            else
                % -------------------- Evaluate Arrhenius factor --------------------
                % Compute Arrhenius factor
                arrhenius_factor = arrhenius_temperature_adjustment(derivative_energy, ...
                                                                    temperature_ref, ...
                                                                    temperature, ...
                                                                    concentration_elyte);
            end

            % -------------------- Evaluate diffusion_coefficient_reference(concentration_elyte, temperature) --------------------
            
            % Initialize diffusion coefficient
            if isnumeric(derivative_coeff_ref)
                % Scalar case
                derivative_coefficient = derivative_coeff_ref;
            elseif isa(derivative_coeff_ref, 'function_handle')
                % Function handle case: check number of input arguments
                nargin_inputs = nargin(derivative_coeff_ref);
                if nargin_inputs == 1
                    derivative_coefficient = derivative_coeff_ref(concentration_elyte);
                elseif nargin_inputs == 2
                    derivative_coefficient = derivative_coeff_ref(concentration_elyte, temperature);
                else
                    error('derivative_activity_coefficient_concentraion_reference_elyte expects 1 or 2 inputs, got %d', n_inputs);
                end
            else
                error('derivative_activity_coefficient_concentraion_reference_elyte must be a scalar or function handle');
            end

            
            % -------------------- Effective diffusion (element-wise) --------------------
            derivative_coefficient_elyte = derivative_coefficient * arrhenius_factor;   

        end

        %% ----- Thermodynamic Factor ----------------------- %%
        function thermodynamic_factor = compute_thermodynamic_factor(obj, ThermalParams, ...
            TransportParams, ...
            concentration_elyte, ...
            temperature)
            
            if nargin < 4 || isnan(temperature) || isempty(temperature)
                temperature =  ThermalParams.cell.('temperature_reference_cell');
            end

            transference_ion_number_elyte = TransportParams.electrolyte.('transference_ion_number_elyte');
            
            effective_conductivity_elyte = obj.compute_effective_conductivity_electrolyte( ...
                ThermalParams, ...
                TransportParams, ...
                concentration_elyte, ...
                temperature);

            derivative_coefficient = TransportParams.electrolyte.('derivative_activity_coefficient_concentraion_reference_elyte');

            if isa(derivative_coefficient, 'function_handle')
                derivative_coefficient = obj.compute_electrolyte_derivative_activity_coefficient( ...
                    ThermalParams, ...
                    TransportParams, ...
                    concentration_elyte, ...
                    temperature);

            elseif isnan(derivative_coefficient)
                derivative_coefficient = 0;
            end

            constants = ConstantParameters();

            thermodynamic_factor = 2 * constants.R_gas * temperature * ...
                               effective_conductivity_elyte * (transference_ion_number_elyte - 1) / ...
                               constants.F .* (1 + derivative_coefficient);

        end
    end
end