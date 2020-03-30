function parsed_inputs = getInputs(inputs,parsed_inputs,input_types)

    inputNames = fieldnames(parsed_inputs);
    N_inputs = length(inputNames);
    
    if nargin==3 && iscell(input_types) % check if permissable data types are provided
        for i = 1:N_inputs
            j = find(strcmp(inputNames{i},inputs));
            if ~isempty(j)
                this_input = inputs{j+1};
                % Check "OR" first:
                useOr = 0;
                for k = 1:size(input_types{i},1)
                    switch input_types{i}{k,1}
                        case 'char'
                            if ischar(this_input)
                                useOr = k;
                                break;
                            end
                        case 'string'
                            if ischar(this_input)
                                useOr = k;
                                break;
                            end
                        case 'struct'
                            if isstruct(this_input)
                                useOr = k;
                                break;
                            end
                        case 'logical'
                            if islogical(this_input)
                                useOr = k;
                                break;
                            end
                        case 'numeric'
                            if isnumeric(this_input)
                                useOr = k;
                                break;
                            end
                        case 'cell'
                            if iscell(this_input)
                                useOr = k;
                                break;
                            end
                        case 'axes'
                            if isgraphics(this_input,'axes')
                                useOr = k;
                                break;
                            end
                        case 'file'
                            if exist(this_input,'file')==2
                                useOr = k;
                                break;
                            end
                        case 'vector'
                            if isvector(this_input)
                                useOr = k;
                                break;
                            end
                    end
                end
                % Next, check "AND":
                failedAnd = 0;
                for k = 1:size(input_types{i},2)
                    switch input_types{i}{useOr,k}
                        case 'char'
                            if ~ischar(this_input)
                                failedAnd = k;
                                break;
                            end
                        case 'struct'
                            if ~isstruct(this_input)
                                failedAnd = k;
                                break;
                            end
                        case 'logical'
                            if ~islogical(this_input)
                                failedAnd = k;
                                break;
                            end
                        case 'numeric'
                            if ~isnumeric(this_input)
                                failedAnd = k;
                                break;
                            end
                        case 'cell'
                            if ~iscell(this_input)
                                failedAnd = k;
                                break;
                            end
                        case 'axes'
                            if ~isgraphics(this_input,'axes')
                                failedAnd = k;
                                break;
                            end
                        case 'file'
                            if ~exist(this_input,'file')==2
                                failedAnd = k;
                                break;
                            end
                        case ''
                            break;
                    end
                end
                if failedAnd
                    error(['Invalid input for ''',inputs{j},'''.'])
                else
                    parsed_inputs.(inputNames{i}) = inputs{j+1};
                end
            end
        end
    else % ignore input data types
        for i = 1:N_inputs
            j = find(strcmp(inputNames{i},inputs));
            parsed_inputs.(inputNames{i}) = inputs{j+1};
        end
    end
    
end