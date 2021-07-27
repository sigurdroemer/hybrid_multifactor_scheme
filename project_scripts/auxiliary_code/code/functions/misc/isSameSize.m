function chk = isSameSize(varargin)
% 
%	Checks if a set of variables are of the same size.
%
% ----------------------------------------------------------------------
% Parameters
% ---------------------------------------------------------------------- 
%   varargin:   [NxM <type>]    Inputs to compare.
%
% ----------------------------------------------------------------------
% Output
% ----------------------------------------------------------------------
%   chk:        [1x1 logical]   True if all are of the same size.
% ----------------------------------------------------------------------                                           

chk = true;
for i=1:length(varargin)-1
    if size(varargin{i},1) ~= size(varargin{i+1},1) ...
            || size(varargin{i},2) ~= size(varargin{i+1},2)
        chk = false;
        return;
    end
end

end