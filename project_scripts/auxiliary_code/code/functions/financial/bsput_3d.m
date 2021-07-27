function p = bsput_3d(s0,K,r,T,sig2,q)
%
% 	Implements the Black-Scholes put price formula but allows for some particular input 
%   dimensions (see below).
%
% -----------------------------------------------------------------------------------------------
% Parameters
% -----------------------------------------------------------------------------------------------
%   s0:     [MxN real]        Asset price.
%   K:      [1x1xL real]      Strike price.
%   r:      [1x1 real]        Interest rate.
%   T:      [1x1 real]        Time-to-maturity.
%   sig2:   [MxN real]        Black-Scholes volatility squared.
%   q:      [1x1 real ...     Dividend yield (default is 0).
%           (optional)] 
%
% -----------------------------------------------------------------------------------------------
% Output
% -----------------------------------------------------------------------------------------------
%   p:      [MxNxL real]    Prices.
% -----------------------------------------------------------------------------------------------

    if ~exist('q','var') || isempty(q); q = 0; end
    c = bscall_3d(s0,K,r,T,sig2,q);
    p = c + K*exp(-r*T) - s0*exp(-q*T);
    
end