classdef Battery < handle
    % Handle class: all modifications happen in-place (no copying)

    properties
        GeometricParams
        ConcentrationParams
        KineticParams
        ThermodynamicParams
        TransportParams
        ElectricalParams
        ThermalParams
        SimulationParams
        AgingParams
        preComputedMatrix
    end

    methods
        function self = Battery(src)
            % Construct from an existing struct (from readExcel) or empty
            if nargin >= 1 && ~isempty(src)
                f = fieldnames(src);
                for k = 1:numel(f)
                    self.(f{k}) = src.(f{k});
                end
            end
        end
        function s = toStruct(self)
            % Snapshot of public properties (fast; value types copy-on-write)
            props = properties(self);
            for k = 1:numel(props)
                s.(props{k}) = self.(props{k});
            end
        end
    end
end