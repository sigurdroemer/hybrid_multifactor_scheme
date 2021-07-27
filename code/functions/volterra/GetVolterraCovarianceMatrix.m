function SIGMA = GetVolterraCovarianceMatrix(K,kappa,dt)
% 
%   Computes the covariance matrix for a set of Gaussian stochastic Volterra integrals all with 
%   respect to the same underlying Brownian motion W(t). Specifically, we return the covariance 
%   matrix of the following random vector Y:
%
%   Y :=   (
%
%                W(dt) - W(0),
%
%                int_{0}^{dt} K_1(dt-s) dW(s),
%                int_{0}^{dt} K_1(2*dt-s) dW(s),
%                ...,
%                int_{0}^{dt} K_1(kappa(1)*dt-s) dW(s)
%                                   
%                ...,
%
%                int_{0}^{dt} K_m(dt-s) dW(s),
%                int_{0}^{dt} K_m(2*dt-s) dW(s),
%                ...,
%                int_{0}^{dt} K_m(kappa(m)*dt-s) dW(s)
%
%                 )'
%
%   Here dt > 0 is the step size and K_1,...,K_m kernel functions.
%
%   Covariances are computed by numerical integration.
%
% ------------------------------------------------------------------------------------------------
%  Parameters
% ------------------------------------------------------------------------------------------------
%   K:      [1xm cell]      Kernel functions. Each element is an object of the type 
%                           [1x1 KernelFunctionClass] representing K_i for i=1,...,m.
%   kappa:  [1xm integer]   Number of shifted Volterra integrals per kernel. 
%   dt:     [1x1 integer]   Step size.
%
% ------------------------------------------------------------------------------------------------
%  Output
% ------------------------------------------------------------------------------------------------
%   SIGMA:  [sum(kappa) + 1 x sum(kappa) + 1 real]  Covariance matrix.
%
% ------------------------------------------------------------------------------------------------

if size(kappa,2) ~= size(K,2) || size(kappa,1) > 1 || size(K,1) > 1
    error('GetVolterraCovarianceMatrix: Invalid input sizes.');
end

% Remove terms where kappa = 0:
idxKeep = kappa ~= 0;
kappa = kappa(idxKeep);
K = K(idxKeep);

% Initialize:
m = size(K,2);
SIGMA = NaN(sum(kappa)+1,sum(kappa)+1);
SIGMA(1,1) = dt;

if sum(kappa) == 0
    return;
end

% Kernel for Brownian part:
Kid = KernelFunctionClass(1,0,@(obj,t)(1),0);

% Compute covariances between Brownian increment and stochastic integrals:
ii = 1;
for i=1:m
    for j=1:kappa(i)
        ii = ii + 1;
        SIGMA(ii,1) = GetVolterraCovariance(Kid,0,dt,dt,K{i},0,dt,j*dt);
    end
end

% Create auxiliary loop vectors:
[idxK,kappa_expanded] = deal(NaN(sum(kappa),1));
ii = 1;
for i=1:size(K,2)
    for j=1:kappa(i)
        idxK(ii) = i;
        kappa_expanded(ii) = j;
        ii = ii + 1;
    end
end

% Compute covariances between stochastic integrals:
for i=1:size(idxK,1)
    for j=1:size(idxK,1)
        if j <= i
            SIGMA(i+1,j+1) = GetVolterraCovariance(K{idxK(i)},0,dt,kappa_expanded(i)*dt,...
                                                   K{idxK(j)},0,dt,kappa_expanded(j)*dt);
        end
    end
end

% Remaining entries follow by symmetry:
SIGMA(isnan(SIGMA)) = 0;
SIGMADiagZero =  SIGMA - diag(diag(SIGMA));
SIGMA = SIGMA + SIGMADiagZero';

end

