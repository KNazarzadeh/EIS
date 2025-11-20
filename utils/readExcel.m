function battery = readExcel(filePath)
    % Load default parameters list
    defaultParametersList = DefaultParameters.defineDefaultParameters();

    % Define the sheets to process
    sheets = sheetnames(filePath);

    % Remove the unwanted sheet named "x"
    sheets = setdiff(sheets, "Experiment Data", 'stable');

    % Initialize the main structure
    battery = struct();
    
    % Loop through each sheet
    for i = 1:length(sheets)
        sheet_name = sheets{i};
        % Create valid field name for the sheet (e.g., GeometricParams)
        sheet_field = [sheet_name, 'Params'];

        % Read the sheet
        try
            if ismember(sheet_name, ["Transport", "Kinetic", "Thermodynamic", "Geometric"])
                opts = detectImportOptions(filePath, 'Sheet', sheet_name);
                opts = setvartype(opts, {'Values'}, 'char');
                data = readtable(filePath, opts, 'Sheet', sheet_name);
            else
                data = readtable(filePath, 'Sheet', sheet_name);
            end
            % data = readcell(filePath, 'Sheet', sheet_name, 'TextType', 'char');
        catch e
            warning('Failed to read sheet %s: %s', sheet_name, e.message);
            continue;
        end
               
        % Get unique domains
        % Remove empty cells or empty strings first
        nonEmptyDomains = data.Domain(~cellfun(@isempty, data.Domain));

        % Then get unique (preserving order)
        domain = unique(nonEmptyDomains, 'stable');

        % domain = unique(data.Domain, 'stable');  % Remove empty strings
        domains = lower(domain(1:end));

        % Loop through each domain
        startIdx = zeros(length(domains),1);
        endIdx = zeros(length(domains),1);

        for i = 1:length(domains)
            domain = domains{i};

            if i == length(domains)
                 idx = find(strcmpi(data.Domain, domain)) + 1;

                 startIdx(i) = idx(1);       % first occurrence index
                 endIdx(i) = length(data.Domain) + 1;
            else
                nextdomain = domains{i+1};
                idx = find(strcmpi(data.Domain, domain)) + 1;
                idxnext = find(strcmpi(data.Domain, nextdomain)) + 1;

                startIdx(i) = idx(1);       % first occurrence index
                endIdx(i) = idxnext - 1;       % last occurrence index
            end
        end

        for i = 1:length(domains)
            domain = domains{i};
            % Determine the target field in the structure
            if strcmp(domain, 'electrode negative')
                target_field = 'electrode.negative';
                % Split the target_field into parts if it has dots (nested structure)
                parts = strsplit(target_field, '.');
                prefix = 'neg';                
            elseif strcmp(domain, 'electrode positive')
                target_field = 'electrode.positive';
                % Split the target_field into parts if it has dots (nested structure)
                parts = strsplit(target_field, '.');
                prefix = 'pos';
            else
                target_field = domain;
                if strcmp(domain, "cell")
                    prefix = "cell";

                elseif strcmp(domain, "separator")
                    prefix = "sep";

                elseif strcmp(domain, "electrolyte")
                    prefix = "elyte";
                end
            end
            % Create a temporary struct for parameters
            paramStruct = struct();
            for idx = startIdx(i):endIdx(i)-1
                % Get parameter names and values for this domain
                paramName = strcat(lower(data.Parameters{idx}), '_', prefix);
                                    
                if contains(paramName, '%')
                    % Replace '%' with ''
                    paramName = strrep(paramName, '%', '');
                end

                if contains(paramName, 'number of')
                    % Replace 'number of ' with 'number'
                    paramName = strrep(paramName, 'number of', 'number');
                end

                % Check if there is a space in the string
                if contains(paramName, ' ')
                    paramName = strrep(paramName, ' ', '_');  % replace space with '_'
                end

                if ismember(sheet_name, ["Transport", "Kinetic", "Thermodynamic"])
                    if contains(paramName, 'equation')
                        % Read as string (text form of function or equation)
                        cellval  = data.Values(idx);
                        if isempty(cellval{1}) || ~isnan(str2double(cellval{1}))
                            paramValue = str2double(cellval{1});
                        elseif isnan(str2double(cellval{1}))
                            if contains(cellval{1}, "stoichiometry", 'IgnoreCase', true)
                                paramValue = str2func(['@(stoichiometry) ', cellval{1}]);
                            elseif contains(cellval{1}, "sto")
                                cellval{1} = replace(cellval{1}, "sto", "stoichiometry");
                                paramValue = str2func(['@(stoichiometry) ', cellval{1}]);
                            end
                        end
                    elseif contains(paramName, ["diffusion", "conductivity", "derivative", "energy"])
                        if ismember(domain, ["electrode negative", "electrode positive"])
                            % Convert formula string to an anonymous function of ce
                            cellval  = data.Values(idx);
                            if isempty(cellval{1}) || ~isnan(str2double(cellval{1}))
                                paramValue = str2double(cellval{1});
                            elseif isnan(str2double(cellval{1}))
                                if contains(cellval{1}, "stoichiometry", 'IgnoreCase', true) && contains(cellval{1}, "T", 'IgnoreCase', true)
                                    paramValue = str2func(['@(stoichiometry, T) ', cellval{1}]);
                                elseif contains(cellval{1}, "stoichiometry", 'IgnoreCase', true)
                                    paramValue = str2func(['@(stoichiometry) ', cellval{1}]);
                                elseif contains(cellval{1}, "sto")
                                    cellval{1} = replace(cellval{1}, "sto", "stoichiometry");
                                    paramValue = str2func(['@(stoichiometry) ', cellval{1}]);
                                end
                            end
                        elseif strcmpi(domain, "electrolyte")
                            % Convert formula string to an anonymous function of ce
                            cellval  = data.Values(idx);
                            if isempty(cellval{1}) || ~isnan(str2double(cellval{1}))
                                paramValue = str2double(cellval{1});
                            elseif isnan(str2double(cellval{1}))
                                paramValue = cellval{1};
                                % paramValue = strrep(paramValue, '/', './'); 
                                % paramValue = strrep(paramValue, '^', '.^');
                                if contains(paramValue, "ce", 'IgnoreCase', true) && contains(cellval{1}, "T", 'IgnoreCase', true)
                                    paramValue = str2func(['@(ce, T) ', paramValue]);
                                else
                                    paramValue = str2func(['@(ce) ', paramValue]);
                                end

                            end
                        end
                    else
                        cellval  = data.Values(idx);
                        paramValue = str2double(cellval{1});
                    end

                elseif ismember(sheet_name, "Geometric")
                    if ismember(domain, ["electrode negative", "electrode positive"])
                        if contains(paramName, ['particle_radius_' prefix])
                            % Read as string (text form the Excel cell)
                            cellval = data.Values(idx);
                            if isnumeric(cellval)
                                paramValue = cellval;
                                number_of_particles = length(paramValue);
                            else
                                strval = cellval{1};
                                paramValue = str2num(strval);
                                number_of_particles = length(paramValue);
                            end
                        elseif contains(paramName, ['particle_fraction_' prefix])
                            % Read as string (text form the Excel cell)
                            cellval = data.Values(idx);
                            if isnumeric(cellval)
                                paramValue = cellval;
                            else
                                strval = cellval{1};
                                paramValue = str2num(strval);
                            end
                        elseif contains(paramName, ['radial_nodes_' prefix])
                            % Read as string (text form the Excel cell)
                            cellval = data.Values(idx);
                            if isnumeric(cellval)
                                paramValue = cellval;
                                number_of_nodes = length(paramValue);
                            else
                                strval = cellval{1};
                                paramValue = str2num(strval);
                                number_of_nodes = length(paramValue);
                            end
                        else
                            paramValue = str2double(data.Values(idx));
                        end
                    else
                        paramValue = str2double(data.Values(idx));
                    end
                else
                    paramValue = data.Values(idx);
                end
                
                if ~isa(paramValue, 'function_handle') && any(isnan(paramValue)) && any(contains(paramName, defaultParametersList))
                    paramValue = DefaultParameters.extractDefaultParameters(material, domain, paramName);
                end

                if ismember(sheet_name, "Thermal")
                    if contains(paramName, "temperature_")
                        paramValue = unitCorrection(paramName, paramValue);
                    end
                end
                % Dynamically add each parameter to the temporary struct
                if ismember(sheet_name, "Geometric")
                    if ismember(domain, ["electrode negative", "electrode positive"])
                        if contains(paramName, ['particle_radius_' prefix])
                            if isvector(paramValue) && ~isscalar(paramValue)
                                for n = 1:number_of_particles
                                    p_name = insertAfter(paramName, "particle_", num2str(n) + "_");
                                    % Initialize paramStruct.particles if needed
                                    if ~isfield(paramStruct, 'particles')
                                        paramStruct.particles = struct();
                                    end
                                    paramStruct.particles.(p_name) = paramValue(n);
                                end
                                paramName = (['number_particles_' prefix]);
                                paramStruct.(paramName) = number_of_particles;
                            else
                                p_name = paramName;
                                 % Initialize paramStruct.particles if needed
                                if ~isfield(paramStruct, 'particles')
                                    paramStruct.particles = struct();
                                end
                                paramStruct.particles.(p_name) = paramValue;
                                paramName = (['number_particles_' prefix]);
                                paramStruct.(paramName) = number_of_particles;
                            end
                        elseif contains(paramName, ['particle_fraction_' prefix])
                            if isvector(paramValue) && ~isscalar(paramValue)
                                for n = 1:number_of_particles
                                    p_name = insertAfter(paramName, "particle_", num2str(n) + "_");
                                    % Initialize paramStruct.particles if needed
                                    if ~isfield(paramStruct, 'particles')
                                        paramStruct.particles = struct();
                                    end
                                    paramStruct.particles.(p_name) = paramValue(n);
                                end
                            else
                                p_name = paramName;
                                 % Initialize paramStruct.particles if needed
                                    if ~isfield(paramStruct, 'particles')
                                        paramStruct.particles = struct();
                                    end
                                paramStruct.particles.(p_name) = paramValue;
                            end
                        elseif contains(paramName, ['radial_nodes_' prefix])
                            if isvector(paramValue) && ~isscalar(paramValue)
                                paramName = (['number_radial_nodes_' prefix]);
                                for n = 1:number_of_nodes
                                    p_name = insertAfter(paramName, "particle_", num2str(n) + "_");
                                    paramStruct.particles.(p_name) = paramValue(n);
                                end
                            else
                                paramName = (['number_radial_nodes_' prefix]);
                                paramStruct.particles.(paramName) = paramValue;
                            end
                        elseif contains(paramName, ['spatial_nodes_' prefix])
                            paramName = (['number_spatial_nodes_' prefix]);
                            paramStruct.(paramName) = paramValue;

                        else
                            paramStruct.(paramName) = paramValue;
                        end
                    else
                        paramStruct.(paramName) = paramValue;
                    end
                else
                    paramStruct.(paramName) = paramValue;
                end
            end
            
            if contains(target_field, '.')
                % Split the target_field into parts if it has dots (nested structure)
                parts = strsplit(target_field, '.');
                % Assign the parameter struct to the target field in the battery struct
                battery.(sheet_field).(parts{1}).(parts{2}) = paramStruct;
            else
                % Assign the parameter struct to the target field in the battery struct
                battery.(sheet_field).(target_field) = paramStruct;
            end

        end
    end

    % Display success message
    fprintf('Battery parameters successfully created from file: \n"%s" \n', filePath);
    disp("------------------------")
end