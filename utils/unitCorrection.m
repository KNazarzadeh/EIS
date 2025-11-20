function param_value = unitCorrection(param_name, param_values)
    if contains(param_name, "temperature")
       if  param_values < 150

           param_value = param_values + 273.15;    % Change Celsius to Klvin [K]
                                                   % Celsius to Kelvin conversion equation: T(K)=T(°C)+273.15
       else
           param_value = param_values;
       end
    end
end
