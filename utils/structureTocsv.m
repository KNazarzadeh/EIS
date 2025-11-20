clear all
clear all
close all
clc

% load 500cycles with STO controller
% outDir = 'C:\Users\k.nazarzadeh\Documents\';
% outFile = fullfile(outDir, 'AgingParams_500_cycles_no_battery.mat');
% AgingParams = load(outFile, 'AgingParams_500_cycles_no_battery');

% fileName = "AgingParams_500_cycles_no_battery.mat";

% AgingParams = AgingParams.AgingParams_500_cycles_no_battery;

% load 500cycles without STO controller
outDir = 'C:\Users\k.nazarzadeh\Documents\';
outFile = fullfile(outDir, 'AgingParams_500cycles_no_sto_control.mat');
AginParams = load(outFile, 'AgingParams_500cycles_no_sto_control');

AgingParams = AgingParams.AgingParams_500cycles_no_sto_control;

fileName = "AgingParams_500cycles_without_sto.mat";


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get the field names
% fieldNames = fieldnames(AgingParams);
% numFields = numel(fieldNames);
% 
% % cycle_num.
% cycle_num = 500;
% 
% % Initialize a new structure to store tables
% CycleData = struct();
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for cycle = 1:cycle_num
%     dataTable = table();
%     % Fill the table with field values
%     for nf = 1:numFields
%         val = AgingParams.(fieldNames{nf}){cycle};
%         dataTable.(fieldNames{nf}) = val(:);
%     end
%     % Keep the original names as a reference (optional)
%     dataTable.Properties.VariableDescriptions = fieldNames;
% 
%     % Store the table in the new structure with field name cycle_1, cycle_2, etc.
%     cycleFieldName = sprintf('cycle_%d', cycle);
%     CycleData.(cycleFieldName) = dataTable;    
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% save(fileName, 'CycleData');

fieldNames = fieldnames(AgingParams);
numFields  = numel(fieldNames);
validNames = matlab.lang.makeValidName(fieldNames,'ReplacementStyle','delete');
cycle_num  = numel(AgingParams.(fieldNames{1}));

CycleData = struct();

for c = 1:cycle_num
    % collect columns and lengths
    cols = cell(numFields,1);
    lens = zeros(numFields,1);
    for nf = 1:numFields
        v = AgingParams.(fieldNames{nf}){c}(:);
        cols{nf} = v;
        lens(nf) = numel(v);
    end
    H = max(lens);              % table height for this cycle

    % build table with padding
    T = table();
    for nf = 1:numFields
        v = cols{nf};
        if numel(v) < H
            if isnumeric(v)
                v(end+1:H,1) = NaN;
            else
                v = string(v); v(end+1:H,1) = missing;
            end
        end
        T.(validNames{nf}) = v;
    end
    T.Properties.VariableDescriptions = fieldNames;  % keep originals

    CycleData.(sprintf('cycle_%d', c)) = T;
end

% save
save(fileName, 'CycleData');
