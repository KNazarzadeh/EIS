classdef ConstantParameters
%
% Standard constants
%
    properties
        R_gas = 8.314;          % Universal gas constant [J/(mol.K)]
        F = 96487;              % Faraday's constant [C/mol]
        T_ref = 25 + 273.15;    % Reference temperature [K]
                                % Celsius to Kelvin conversion equation: T(K)=T(°C)+273.15
        reduced_planck = 1.055E-34; % Reduced Planck constant [J.s]

    end
end
  

