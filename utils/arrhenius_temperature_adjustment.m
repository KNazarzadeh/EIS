function factor = arrhenius_temperature_adjustment(Ea, T_ref, temperature, concentration_elyte)
        % Compute the Arrhenius temperature dependence factor.
        % Inputs:
        % -------
        %   Ea: Activation energy (J/mol)
        %   T_ref: Cell reference temperature (K)
        %   temperature: Temperature (K)

        % Output:
        % -------
        %   factor: Dimensionless
        
        % Load standard constants for
        % R_gas: Universal gas constant ()

        constants = ConstantParameters();

        if isempty(temperature)
            temperature = T_ref;
        elseif temperature <= 0
            error('Temperature must be positive (K).');
        end
        
        if isnumeric(Ea) 
            factor = exp(Ea / constants.R_gas * (1 / T_ref - 1 / temperature));
        elseif isa(Ea, 'function_handle')
            n = nargin(Ea);
            if n == 1
                factor = @(concentration_elyte) exp(Ea(concentration_elyte) / constants.R_gas * (1 / T_ref - 1 / temperature));
                factor = factor(concentration_elyte);
            elseif n == 2
                factor = @(concentration_elyte, temperature) exp(Ea(concentration_elyte, temperature) / constants.R_gas * (1 / T_ref - 1 / temperature));
                factor = factor(concentration_elyte, temperature);
            end

            % paramName = regexp(func2str(Ea), '@\((.*?)\)', 'tokens');
            % paramName = paramName{1}{1};
            % % Validate paramName
            % if ~isvarname(paramName)
            %     error('Invalid parameter name: %s', paramName);
            % end
            % % Define factor dynamically using eval
            % eval(sprintf('factor = @(%s, T) exp(Ea(%s) ./ %f * (1 / 298 - 1 ./ T));', paramName, paramName, constants.R_gas));
        end

end