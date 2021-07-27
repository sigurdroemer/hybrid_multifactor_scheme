function ok = VerifyStructFieldDimensions(X,fieldNames,expectedDims,allowEmpty,allowNonExistent)
%   
%   Checks if fields of a given struct have the expected dimensions.
%
% ------------------------------------------------------------------------------------------------
%  Parameters
% ------------------------------------------------------------------------------------------------
%  X:                [1x1 struct]    Struct to verify fields for.
%  fieldNames:       [1xN cell]      Field names to check.
%  expectedDims:     [1xN cell]      Expected dimensions of corresponding field. Each element
%                                    must either be a string with options 'row_vector',
%                                    'column_vector' or a [1x2 integer] vector specifying the 
%                                    expected dimensions of the matrix.
%  allowEmpty:       [1xN logical]   If true we allow the field to be empty.
%  allowNonExistent: [1xN logical]   If true we allow the field to not exist. I.e. we only 
%                                    check the dimensions if the field exists.
%
% ------------------------------------------------------------------------------------------------
%  Output
% ------------------------------------------------------------------------------------------------
%	ok:		[1x1 logical]	True if fields are valid, otherwise false.
%
% ------------------------------------------------------------------------------------------------

ok = true;
for i=1:size(fieldNames,2)
    % Check if field even exists:
    if ~allowNonExistent(i) && ~isfield(X,fieldNames{i})
        ok = false;
        return;
    elseif ~isfield(X,fieldNames{i})
        continue;
    end
    % Check if field is empty:
    if ~allowEmpty(i) && isFieldEmpty(X,fieldNames{i})
        ok = false;
        return;
    elseif isFieldEmpty(X,fieldNames{i})
        continue;
    end
    % Check dimensions of field:
    dimActual = size(X.(fieldNames{i}));    
    if isstring(expectedDims{i}) || ischar(expectedDims{i})
        if strcmpi(expectedDims{i},'row_vector') && dimActual(1) > 1
            ok = false;
            return;
        elseif strcmpi(expectedDims{i},'column_vector') && dimActual(2) > 1
            ok = false;
            return;
        end
    else
        if max(abs(dimActual-expectedDims{i})) > 0
            ok = false;
            return;
        end
    end
    
end

end

