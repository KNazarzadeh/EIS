classdef DataPreprocessing
    % DataPreprocessing: Preprocess data for the simulation

    methods(Static)
        % function obj = DataPreprocessing(timeRaw, input_current, timeStep)
        %     % Constructor: Preprocess time and current, compute all outputs
        %     % obj.timeRaw = timeRaw;
        %     % obj.input_time = obj.preprocess_RawTime(timeRaw);
        %     % obj.timeVector = obj.define_timeVector(obj.input_time, timeStep);
        %     % obj.applied_current = obj.compute_applied_current(input_current, obj.input_time, obj.timeVector);
        % end 

        function applied_current = preprocessCurrent(battery, applied_current, isCCCV, ...
                dischargeC_rate, chargeC_rate, cycle_startmode, parallelCellNum)
            timeFinal = battery.SimulationParams.timeFinal;
            if isCCCV
                if strcmpi(cycle_startmode, "discharge")
                    applied_current = -dischargeC_rate * battery.ElectricalParams.cell.('nominal_capacity_cell');
                elseif strcmpi(cycle_startmode, "charge")
                    applied_current = chargeC_rate * battery.ElectricalParams.cell.('nominal_capacity_cell');
                end
                applied_current = applied_current/parallelCellNum;
            
            else
                if ~isvector(applied_current)
                    if ~isfield(battery.ElectricalParams.cell, 'nominal_capacity_cell')
                        error('Nominal capacity required for C-rate calculation.');
                    end        
                elseif isnan(applied_current)
                    if ~isfield(battery.ElectricalParams.cell, 'nominal_capacity_cell')
                        error('Nominal capacity required for C-rate calculation.');
                    end

                elseif isvector(applied_current) && ~isscalar(applied_current)
                    applied_current = applied_current;

                elseif isscalar(applied_current)
                    if strcmpi(cycle_startmode, "charge") && applied_current > 0
                        applied_current = applied_current * ones(timeFinal,1);
                    elseif strcmpi(cycle_startmode, "discharge") && applied_current < 0
                        applied_current = applied_current * ones(timeFinal,1);
                    end
                end
            end
        end

        function input_time = preprocess_RawTime(obj, timeRaw)
            % Preprocess time: Normalize to start at 0, ensure monotonicity
            if isempty(timeRaw)
               error('Experiment time is empty!');
            end
            
            if timeRaw(1) ~= 0
                % normalize time to start from 0 second
                input_time = timeRaw - timeRaw(1);
            else
                input_time = timeRaw;
            end
            
            if any(diff(input_time) <= 0)
                error('Time vector must be strictly increasing.');
            end
        end

        function timeVector = define_timeVector(obj, input_time, timeStep)
            % Create time vector with specified time step
            if isempty(input_time)
               error('Time is empty!');
            end
            timeVector = (0:timeStep:input_time(end))';
        end

        function applied_current = compute_applied_current(obj, input_current, input_time, timeVector)
            % Interpolate current at timeVector
            if isempty(timeVector)
               error('Time is empty!');
            end
            
            applied_current = interp1(input_time, input_current, timeVector);

        end
    end
end