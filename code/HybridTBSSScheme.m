function [Y,t,dW] = HybridTBSSScheme(N,n,T,sigma,alpha,kappa,Z,conv_method,W_bold,cov_method,SIGMA)
% Description: Simulates the truncated Brownian semistationary (TBSS) process
%
%   Y(t) = int_0^t g(t-s)*sigma dW(s)
%
% where sigma is a constant, W(t) a Brownian motion and g(x) = x^(alpha) with -0.5 < alpha < 0.5, 
% alpha <> 0.
%
% The function uses the hybrid scheme from (Bennedsen et al., 2017).
%
% Remarks:
%   o Some input parameters require knowing floor(n*T). In the code this is computed as 
%           floor(n*T+eps(n*T))                                                             (1)
%     and not as
%           floor(n*T).                                                                     (2)
%     This is to avoid round off errors. Whenever floor(n*T) is written here in the description
%     it is exactly the numerical computation (1) that is meant.
%   o If the covariance matrix used for the hybrid scheme is computed to not be positive definite 
%     (can happen e.g. due to round off errors) we use the modified Cholesky factorisation   
%     algorithm from (Cheng and Higham, 1998) with the code obtained from 
%     https://github.com/higham/modified-cholesky (retrieved on the 17th of May 2020).
% 
% Required parameters:
%   N:     [1x1 integer] Number of samples to generate.
%   n:     [1x1 integer] Number of steps per unit of time, i.e. steps are of the length dt := 1/n.
%   T:     [1x1 real] Final time point to simulate till. To fit the discretisation the process 
%           will then be simulated on the M := floor(n*T) + 1 grid points 
%           0 < 1/n < 2/n < ... < floor(n*T)/n                                          (3)
%   sigma: [1x1 real] See the description.
%   alpha: [1x1 real] See the description.
%   kappa: [1x1 integer] An integer between 0 and floor(n*T) specifying the number of 
%          sub-integrals that should be simulated exactly. Corresponds to the kappa-parameter 
%          from (Bennedsen et al., 2017).
%
% Optional parameters (leave empty as [] if not to be used):
%   Z:            [(M-1)*N x (kappa+1) real or empty] i.i.d. standard normals needed for the 
%                 simulation. Can be left unspecified or empty in which case they are automatically 
%                 simulated in the code.
%
%   conv_method:  [1x1 string or empty] Options are 'fft' and 'conv_loop'. If 'fft' is 
%                 chosen we use the fast Fourier transform to compute the convolution. The choice 
%                 'conv_loop' computes the convolution via an explicit loop in the time dimension 
%                 (not recommended). The default is 'fft'.
%
%   W_bold:       [N x M - 1 x kappa + 1 real or empty] Random variables (Gaussian) to be 
%                 used by the scheme after transforming. Cannot be used at the same time as the 
%                 'Z' parameter. Should be used with care to ensure paths are simulated correctly. 
%                 Inspect the code for how this matrix should be generated. If left empty the 
%                 matrix is automatically computed/sampled in the code.
%
%   cov_method:
%                 Methods for computing the covariance matrix of the scheme. Options are 'exact' 
%                 in which case the matrix is computed exactly using the 'hypergeom' function and 
%                 'numerical_integration' in which case the covariances are computed by numerical 
%                 integration. The former choice can be very slow for large values of kappa. The
%                 default value is therefore 'exact' for kappa <= 5 and otherwise
%                 'numerical_integration'.
%
%   SIGMA:        [(kappa+1) x (kappa+1) real] Covariance matrix of the W_bold random variables.
%                 It will be automatically computed in the code if left unspecified (also the 
%                 default).
%                 
% Output:
%   Y:   [NxM real] Simulated values of Y(t).
%   t:   [1xM real] Time points.
%   dW:  [Nx(M-1) real] Increments of the underlying Brownian motion W(t) across the time points
%        from (3).
%
% References:
%   o Bennedsen, M., Lunde, A. and Pakkanen, M.S., Hybrid scheme for Brownian semistationary 
%     procesess. Finance and Stochastics, 2017, 21(4), 931-965.
%   o Cheng, S.H., and Higham, N.J., A modified Cholesky algorithm based on a symmetric indefinite 
%     factorization. SIAM Journal on Matrix Analysis and Applications, 1998, 19(4), 1097-1110.
%

% The below computation avoids the round off errors that can arise using floor(n*T):
floor_nT = floor(n*T+eps(n*T));
t = (0:1:floor_nT)'/n;
M = floor_nT + 1;

if ~isempty(Z) && ~isempty(W_bold)
    error(['HybridTBSSScheme: You cannot use both the ''Z'' and ''W_bold'' parameter at the ',...
           'same time.']);
end

% Check alpha value is ok:
if alpha < -0.5 || alpha > 0.5 || alpha == 0
    error('HybridTBSSScheme: ''alpha'' parameter must be between -0.5 and 0.5, zero excluded.');
end

% Check kappa value is ok:
if mod(kappa,1) ~= 0
    error('HybridTBSSScheme: Kappa must be integer valued.');
end
if kappa < 0 || kappa > floor_nT
    error('HybridTBSSScheme: Kappa must be an integer between 0 and floor(n*T).');
end
    
if isempty(W_bold)
    % Simulate random variables:
    if isempty(Z)
        Z = randn((M-1)*N,kappa+1);
    else
        if size(Z,1) ~= (M-1)*N || size(Z,2) ~= kappa + 1
            error(['HybridTBSSScheme: The Z-matrix containing the i.i.d. standard normal ',...
                   'random variables must be of the size [(M-1)*N x kappa + 1].']);
        end
    end
    
    % Compute covariance matrix:
    if ~exist('covmat_method','var') || isempty(cov_method)
        if kappa <= 5
            cov_method = 'exact';
        else
            cov_method = 'numerical_integration';
        end
    end
    if exist('SIGMA','var') && ~isempty(SIGMA)
        if size(SIGMA,1) ~= kappa + 1 || size(SIGMA,2) ~= kappa + 1
            error(['HybridTBSSScheme: Invalid dimensions of ''SIGMA'' matrix. Dimensions ',...
                 'must be (kappa+1) x (kappa+1).']);
        end
    else
        if strcmpi(cov_method,'exact')
            SIGMA = CovMatrixHybrid(n,kappa,alpha,'double');
        elseif strcmpi(cov_method,'numerical_integration')
            K = @(t)(t.^alpha);
            SIGMA = CovMatrixByNumericalIntegration(K,kappa,n);
        else
            error('HybridTBSSScheme: Invalid choice of parameter ''covmat_method''.');
        end
    end

    try
        % Attempt a Cholesky factorisation:
        A = chol(SIGMA, 'lower');
    catch
        % If the matrix is not positive definite:
        [L, D, P] = modchol_ldlt(SIGMA);
        SIGMA_pert = P'*L*D*L'*P;

        % Check the pertubation matrix:
        E = SIGMA_pert - SIGMA;
        max_pert = max(max(abs(E)));
        warning(['HybridTBSSScheme: The covariance matrix was not positive ',...
                 'definite. This could be due to numerical round off errors. Using ',...
                 'modified Cholesky factorisation instead. The largest entry of the pertubation'...
                 ' matrix is ', num2str(max_pert),'.']);

        A = chol(SIGMA_pert, 'lower');
    end

    W_bold = reshape(Z*(A.'),N,M-1,kappa+1);
end
dW = squeeze(W_bold(:,:,1));

% Initialize the Y-matrix:
Y = zeros(N,M);
Y(:,1) = 0;

% Compute sums involving the Volterra integrals:
for i=1:M-1
    temp = 0;
    for k=1:min(i,kappa)
        temp = temp + sigma.*W_bold(:,i-k+1,k+1);
    end
    Y(:,i+1) = temp;
end

% Compute weights for the convolution:
ksUpper = (kappa+1:floor_nT);
bk = ((ksUpper.^(alpha + 1) - (ksUpper - 1).^(alpha + 1))/(alpha + 1)).^(1/alpha);    
weights = (bk/n).^(alpha);

if isempty(conv_method)
    conv_method = 'fft';
end

% Compute convolution:
if strcmpi(conv_method,'fft')
    n_elem = size(weights,2);
    n_padded = 2^nextpow2(2*n_elem-1);
    Y2 = ifft(fft(weights,n_padded).*fft(sigma.*dW(:,1:M-1-kappa),n_padded,2),[],2);        
elseif strcmpi(conv_method,'conv_loop')
    Y2 = zeros(N,M-1-kappa);
    weights = fliplr(weights);
    for i=1:size(Y2,2)
        Y2(:,i) = sum(weights(end-i+1:end).*sigma.*dW(:,1:i),2);
    end
end

% Add convolution to the result matrix:
Y(:,kappa+2:end) = Y(:,kappa+2:end) + Y2(:,1:M-1-kappa);

end

