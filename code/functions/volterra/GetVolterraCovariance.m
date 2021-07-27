function covVal = GetVolterraCovariance(K1,a1,b1,c1,K2,a2,b2,c2)
%
%   Computes the covariance between
%
%       int_a1^b1 K1(c1-s) dW(s)         and         int_a2^b2 K2(c2-s) dW(s)
%
%   where W(t) is a Brownian motion.
%
%   We assume 0 <= a1 < b1 <= c1, and similarly for the second stochastic integral.
%
% -----------------------------------------------------------------------------------------------
%  Parameters
% -----------------------------------------------------------------------------------------------
%   K1:         [1x1 KernelFunctionClass]   See description.
%   a1:         [1x1 real]                  See description.
%   b1:         [1x1 real]                  See description.
%   c1:         [1x1 real]                  See description. 
%
%   K2:         [1x1 KernelFunctionClass]   See description.
%   a2:         [1x1 real]                  See description.
%   b2:         [1x1 real]                  See description.
%   c2:         [1x1 real]                  See description.
%
% -----------------------------------------------------------------------------------------------
%  Output
% -----------------------------------------------------------------------------------------------
%   covVal:     [1x1 real]                  Covariance.
%
% -----------------------------------------------------------------------------------------------

if ~( 0 <= a1 && a1 < b1 && b1 <= c1 &&  0 <= a2 && a2 < b2 && b2 <= c2)
    error('GetVolterraCovariance: Invalid inputs.');
end

a = max(a1,a2);
b = min(b1,b2);

if b <= a
    % Intervals do not overlap => zero covariance:
    covVal = 0;
    
else
    % We use the Ito isometry:
    if c1 == b && c2 == b
        % Both kernels are used in their singular parts:
        [~,~,~,K1K2] = ExpandProductOfKernels(K1,K2);   
        covVal = K1K2.IntegrateVolterra(a,b,b,1,@(s)(1));
        
    elseif c1 == b && c2 ~= b
        % Only the first kernel is used around its singular part
        covVal = K1.IntegrateVolterra(a,b,c1,1,@(s)(K2.Eval(c2-s)));
        
    else
        % Either second kernel is used at singular part or both kernels are non-singular:
        covVal = K2.IntegrateVolterra(a,b,c2,1,@(s)(K1.Eval(c1-s)));
        
    end
	
end


end

