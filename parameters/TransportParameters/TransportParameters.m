classdef TransportParameters
    % Extract and organize transport parameters for battery models
    %
    % Description:
    % -----------
    % This class extracts transport parameters from a specified battery type
    % (e.g., NMC811_chen2020) and organizes them into domains (cell, electrode,
    % electrolyte, separator, sei) based on the electrochemical model type
    % (e.g., SPM, ESPM, DFN).
    %
    methods(Static)
        %---------------------------------------------------------------
        function initialTransportParameters(battery)

            % Extract and organize transport parameters by domain
            %
            % Description:
            % ------------
            % Instantiates the battery type class (e.g., NMC811_chen2020) and extracts
            % transport parameters, organizing them into substructures: cell, electrode
            % (negative, positive, including particle parameters), electrolyte, separator,
            % and sei. For DFN, includes node-based parameters for spatially resolved transport,
            % assuming uniform properties across nodes.
            %
            % Outputs:
            % --------
            % transportParams: Updated object with transportParams structure populated

            % --------------------------------------------------------------------
            % 1. Initialise calculator & aging flags
            % --------------------------------------------------------------------
            calculator = TransportParametersCalculator();

            % Initialize transport parameters structure with domains
            
            % Cell domain
            % No transport parameters in cell domain
            
            % Electrode domain
            % No transport parameters in Electrode domain
            
            % --------------------------------------------------------------------
            % Model-specific geometry (executed once, after electrodes)
            % --------------------------------------------------------------------
            % ESPM Model: Extended Single Particle Model (ESPM) uses electrolyte and Bruggeman parameters
            battery.TransportParams.electrolyte.brug_correction_elyte = ...
                calculator.compute_brugCorrection(battery);
        end
    end
end