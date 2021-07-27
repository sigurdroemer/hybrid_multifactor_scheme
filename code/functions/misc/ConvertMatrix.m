function x = ConvertMatrix(x,prec)
% 
%   Converts a numeric matrix to the specified precision.
%
% ------------------------------------------------------------------------------------------------
%   Parameters
% ------------------------------------------------------------------------------------------------
%      x:       [NxM real]      Numeric matrix to be converted.
%   prec:       [1x1 string]    Options are 'double' and 'single'.
%
% ------------------------------------------------------------------------------------------------
%   Output
% ------------------------------------------------------------------------------------------------
%    x:         [NxM real]      Numeric matrix converted.
%
% ------------------------------------------------------------------------------------------------

   if strcmpi(class(x),'double') && strcmpi(prec,'single')
       x = single(x);
   elseif strcmpi(class(x),'single') && strcmpi(prec,'double')
       x = double(x);
   end
   
end