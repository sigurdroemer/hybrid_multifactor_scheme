function chk = isFieldNonEmpty(x,fieldName)
% 
%	Returns true for a struct if the particular field (1) exists and (2)
%	is non-empty. Returns false otherwise.
%
% -------------------------------------------------------------------------------
% 	Parameters
% -------------------------------------------------------------------------------
%	x:			[1x1 struct] Struct.
%	fieldName:	[1x1 string] Field name to look for.
%
% -------------------------------------------------------------------------------
%   Output
% -------------------------------------------------------------------------------
%	chk:	[1x1 logical] See main description.
%
% -------------------------------------------------------------------------------

chk = isfield(x,fieldName) && ~isempty(x.(fieldName));

end

