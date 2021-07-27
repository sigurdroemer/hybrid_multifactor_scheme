function val = IntegrateFractional(beta,f,a,b,delta)
%
% 	Numerically computes an integral of the form 
%
%       int_a^b s^(beta).*f(s) ds                                                       (*)
%
%   where f is non-singular on the interval [a,b].
%
%   A remark on the computational details:
%   Letting delta > 0 be a (very) small positive number and defining c = max(min(delta,b),a), 
%   i.e. c being delta truncated to the interval [a,b], we decompose the integral (*) as
%
%       int_a^b s^(beta).*f(s) ds = int_a^c s^beta f(s) ds + int_c^b s^beta f(s) ds 
%
%                                 =: I_1 + I_2.
%
%   Here I_2 is evaluated using Matlab's integral function. I_1 is instead evaluated analytically
%   as follows:
%   
%       I_1 approx(=) 0.5*(f(a)+f(c))*int_a^c s^beta ds 
%                  =  0.5*(f(a)+f(c))*(c^(beta+1) - a^(beta+1))/(beta+1).
%
%   * Intended use of parameter delta is to avoid using Matlab's 'integral' function near the 
%     potential singularity since that function can be very slow (and even very inaccurate) when 
%     the integrand is singular.
%
% -----------------------------------------------------------------------------------------------  
%  Parameters
% -----------------------------------------------------------------------------------------------  
%   beta:       [1x1 real]          Fractional index.
%   f:          [1x1 function]      A non-singular function. Must be vectorised.
%   a:          [Nx1 or 1xN real]   Lower integration end point(s).
%   b:          [Nx1 or 1xN real]   Upper integration end point(s).
%   delta:      [1x1 real]          A (very) small positive number. See the description.
% -----------------------------------------------------------------------------------------------  
% Output
% -----------------------------------------------------------------------------------------------  
%   val:        [Nx1 or 1xN real]      Value of integral(s).
% -----------------------------------------------------------------------------------------------  

    % Vectorise:
    val = zeros(size(a));
    for i=1:numel(val)
        val(i) = IntegrateFractionalSub(beta,f,a(i),b(i),delta);
    end

end
    
function val = IntegrateFractionalSub(beta,f,a,b,delta)

    c = max(min(delta,b),a);
    val = 0;
    
    if c > a
        val = val + 0.5*(f(a)+f(c))*(c.^(beta+1) - a.^(beta+1))./(beta + 1) ;
    end
    
    if b > c
        val = val + integral(@(s)((s.^beta).*f(s)),c,b);
    end
    
end




