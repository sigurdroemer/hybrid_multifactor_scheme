function idx = LocateTimepointsInDiscretisation(tDiscr,tLocate)
%
%   Given a discretisation of time points (tDiscr), we for each timepoint in tLocate find 
%   the nearest grid point at or below the value.
%
%   Remark: Round off errors are accounted for.
%
% ------------------------------------------------------------------------------------------------
%  Parameters
% ------------------------------------------------------------------------------------------------
%   tDiscr:         [1xN real]      Must be sorted in ascending order. Values must be unique.
%   tLocate:        [1xM real]      Timepoints to locate in tDiscr by truncating down to the  
%                                   nearest value. We require min(tDiscr) <= min(tLocate).
% ------------------------------------------------------------------------------------------------
%  Output
% ------------------------------------------------------------------------------------------------
%   idx: [1xM integer]  Indices for time points in discretisation.
% ------------------------------------------------------------------------------------------------

% Validate inputs:
if size(tDiscr,1) > 1 || size(tLocate,1) > 1
    error('LocateTimepointsInDiscretisation: Inputs must be row vectors.');
end

if size(unique(tDiscr),2) ~= size(tDiscr,2)
    error('LocateTimepointsInDiscretisation: Values in ''tDiscr'' vector must be unique.');
end

if ~issorted(tDiscr,'ascend')
    error(['LocateTimepointsInDiscretisation: Values in ''tDiscr'' vector must be sorted ',...
           'in ascending order.']);
end

% Remark: We use eps(...) to account for round off errors.
if min(tLocate) + eps(min(tLocate)) < min(tDiscr)
    error('LocateTimepointsInDiscretisation: We require min(tDiscr) <= min(tLocate).');
end

% Locate time points:
idx = zeros(size(tLocate));
for i=1:size(idx,2)
    idx(i) = find(tDiscr <= tLocate(i) + eps(tLocate(i)),1,'last'); 
end

end

