function chk = isStructFieldsSameSize(x)
% 
%	Checks if all fields of a struct are of the same size. 
%
% ----------------------------------------------------------------------------------------
% 	Parameters                                           
% ----------------------------------------------------------------------------------------
%   x:     [1x1 struct]    Struct to check. Fields must be of the form [NxM <type>].
%
% ----------------------------------------------------------------------------------------
%   Output                                            
% ----------------------------------------------------------------------------------------
%   chk:   [1x1 logical]   True if all fields are of the same size.
%
% ----------------------------------------------------------------------------------------

chk = true;
fn = fieldnames(x);
for j=1:size(fn,1)
	if j==1
		n = size(x.(fn{j}),1);
		m = size(x.(fn{j}),2);
	else
		if size(x.(fn{j}),1) ~= n || size(x.(fn{j}),2) ~= m
			chk = false;
			return;
		end
	end
end

end

