function SIGMA = CovMatrixByNumericalIntegration(K,kappa,n)
% Description: Returns the covariance matrix needed for the hybrid-exponential scheme when the 
% SVE to be simulated is one-dimensional and uses the same kernel for both the drift term and 
% the Brownian term. Covariances are obtained by numerical integration.
%
% Parameters:
%   K:      [1x1 function] Kernel function. Must be vectorized.
%   kappa:  [1x1 integer] The number of sub-integrals to simulate exactly in each iteration for 
%           the hybrid multifactor scheme. Corresponds to the same parameter from the function 
%           HybridMultifactorScheme.
%   n:      [1x1 integer] The number of steps per unit of time.
%
% Output: 
%   SIGMA:  [(kappa+1)x(kappa+1) real] Covariance matrix needed for the scheme.
%

SIGMA = NaN(kappa+1,kappa+1);
SIGMA(1,1) = 1/n;

if kappa == 0
    return;
end

for j=2:kappa+1
    SIGMA(1,j) = integral(K,(j-2)/n,(j-1)/n);
end

for k=2:kappa+1
    for j=2:k
        SIGMA(j,k) = integral(@(s)(K(s).*K(s+(k-j)/n)),(j-2)/n,(j-1)/n);
    end
end

SIGMA(isnan(SIGMA))=0;
SIGMADiagZero =  SIGMA - diag(diag(SIGMA));
SIGMA = SIGMA + SIGMADiagZero';

end

