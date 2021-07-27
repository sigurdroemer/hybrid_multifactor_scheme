function Lambda = GetVolterraCovarianceMatrix2(K1,K2,rho,a,b,c,idxK)
%
%   Computes covariance matrix of the following random vector
%
%          int_{a(i)}^{b(i)} K_{idxK(i)}(c(i)-s) dW_{idxK(i)}(s),    i=1,...,N.
%
%   Here dW_1(t)dW_2(t) = rho dt for a rho in [-1,1].
%
% ------------------------------------------------------------------------------------------------
%  Parameters
% ------------------------------------------------------------------------------------------------
%   K1:         [1x1 KernelFunctionClass]   First kernel.
%   K2:         [1x1 KernelFunctionClass]   Second kernel.
%   rho:        [1x1 real]                  Correlation between two Brownian motions.
%   a:          [Nx1 real]                  See description.
%   b:          [Nx1 real]                  See description.
%   c:          [Nx1 real]                  See description.
%   idxK:       [Nx1 integer]               Indices of kernels. K1 ~ 1, K2 ~ 2.
%
% ------------------------------------------------------------------------------------------------
%  Output
% ------------------------------------------------------------------------------------------------
%   Lambda:     [NxN real] Covariance matrix.
%
% ------------------------------------------------------------------------------------------------

N = size(a,1);
Lambda = NaN(N,N);
K = {K1,K2};

% Compute covariances between stochastic integrals:
for i=1:size(idxK,1)
    for j=1:size(idxK,1)
        if j <= i
            Lambda(i,j) = GetVolterraCovariance(K{idxK(i)},a(i),b(i),c(i),...
                                                K{idxK(j)},a(j),b(j),c(j)); 
            if idxK(i) ~= idxK(j)
                Lambda(i,j) = rho*Lambda(i,j);
            end
        end
    end
end

% Remaining entries follow by symmetry:
Lambda(isnan(Lambda)) = 0;
SIGMADiagZero =  Lambda - diag(diag(Lambda));
Lambda = Lambda + SIGMADiagZero';


end

