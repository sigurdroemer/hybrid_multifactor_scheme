function I = SimpsonsIntegration(F,a,b,N)
%
%   Implements Simpson's 1/3 integration rule.
%
% ------------------------------------------------------------------------------------------------
%  Parameters
% ------------------------------------------------------------------------------------------------
%   F:  [1x1 function]  Multivariate function.
%   a:  [1x1 real]      Lower interval point.
%   b:  [1x1 real]      Upper interval point.
%   N:  [1x1 integer]   Number of integration points. Must be odd.
%
% ------------------------------------------------------------------------------------------------
%  Output
% ------------------------------------------------------------------------------------------------
%   I:      [nx1 real]      Integral values for each dimension of F:
% ------------------------------------------------------------------------------------------------

if mod(N,2) ~= 1
    error('SimponsIntegration: ''N'' must be odd.');
end

% Setup:
dx = (b - a)/(N-1);
x = (a:dx:b);

% Compute function values:
Fx = F(x);

% Compute integrals:
I = sum((1/3)*dx*(Fx(:,1:2:end-2) + 4*Fx(:,2:2:end-1) + Fx(:,3:2:end)),2);

end

