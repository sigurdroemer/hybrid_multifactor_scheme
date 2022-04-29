function [Xpaths,Upaths,dW] = HybridMultifactorScheme(N,n,T,g,b,sigma,gamm,c,kappa,varargin)
%
%   Simulates the stochastic Volterra equation
%
%       X(t) = g(t) + int_0^t K(t-s)b(s)ds + int_0^t K(t-s)sigma(s)dW(s)
%
%   where b(t) = b(t,X(t)) and sigma(t) = sigma(t,X(t)) are functions (or alternatively, 
%   pre-simulated processes) and W(t) is a Brownian motion.
%
%   Function uses the hybrid multifactor scheme as presented in (Roemer, 2021). Here the K(t)
%   function is used exactly (i.e. without approximation) for a number of 'kappa' intervals near the 
%   origin and approximated by a sum of exponentials on the remainder of the domain. Specifically, 
%   letting 'dt' be the step size, we for t >= kappa*dt approximate
% 
%     K(t) approx(=) sum_{j=1}^m c_j*exp(-gamm_j*t)
%
%   where c = (c_1,...,c_m)' and gamm = (gamm_1,...,gamm_m)' are appropriate coefficients.
%
%   Remarks: 
%    o Function has been optimised for speed. Code may thus be less readable in some places.
%    o Some input parameters require knowing floor(n*T). We here use the function floor2 that
%      accounts for round off errors. Whenever floor(n*T) is written here in the description
%      it is exactly the numerical computation floor2(n*T) that is meant.
% 
% -------------------------------------------------------------------------------------------------
%  Main parameters:
% -------------------------------------------------------------------------------------------------
%   N:     [1x1 integer] Number of samples to generate.
%   n:     [1x1 integer] Number of steps per unit of time, i.e. steps are of the length dt := 1/n.
%   T:     [1x1 real or empty] Final time point to simulate till. The process will then be 
%          simulated on the M := floor(n*T) + 1 grid points stated below:
%
%               0 < 1/n < 2/n < ... < floor(n*T)/n.                                         (3)
% 
%          The parameter can only be left empty if both the 'tX' and 'tU' parameters are instead 
%          specified. See the description under 'Additional parameters' for an explanation of how  
%          to use those parameters instead.
%   g:     [1x1 function, NxM real or 1x1 real] Parameter represents the g(t) function from the 
%          main description. Parameter is interpreted differently depending on the input type:
%
%               o [1x1 function]:  A vectorised function of one variable representing g(t).
%               o [NxM real]:      Pre-computed values of g(t). We here allow different values
%                                  for different simulation paths. Matrix represents values for 
%                                  all N paths and across all the M time points: 
%                                  0,1/n,...,floor(n*T)/n.
%               o [1x1 real]:      We interpret g(t) equal to this constant value.
%
%   b:     [1x1 function, Nx(M-1) real or 1x1 real] Parameter represents the b(t,x) function from 
%          the main description. Parameter is interpreted differently depending on the input type:
%
%               o [1x1 function]:  A function of two variables representing b(t,x). Must be 
%                                  vectorised in the second argument.
%               o [Nx(M-1) real]:  Pre-simulated values of b(t,x) for all N paths and across 
%                                  all the M - 1 time points: 0,1/n,...,(floor(n*T)-1)/n.
%               o [1x1 real]:      We interpret b(t,x) equal to this constant value. 
%                                  Recommended if you want b(t,x) constant as code has been 
%                                  further optimised for this case.
%
%          Remark: While the solution of your particular SVE may be positive in theory, this
%          will not automatically hold for the numerical scheme. If you specify 'b' as a 
%          function you should therefore either (1) check that the function returns a valid 
%          value even if the x input is negative or (2) use the 'positive' parameter described 
%          further down.
%
%   sigma: [1x1 function, Nx(M-1) real or 1x1 real] As for the 'b' parameter but for sigma(t,x).
%          See also the remark about positivity under 'b'.
%   gamm:  [mx1 real] Exponents for the exponential approximation. See the main description.
%   c:     [mx1 real] Coefficients for the exponential approximation. See the main description.
%   kappa: [1x1 integer] Number of sub-integrals to approximate K exactly. Thus kappa = 0 
%          corresponds to a pure exponential approximation of the kernel function.
%
% -------------------------------------------------------------------------------------------------
%  Additional parameters: 
% -------------------------------------------------------------------------------------------------
% Must be given in name-value pairs. Also, all parameters below are optional except for a few
% special cases to be explained:
%
%   o Either the 'K' parameter must be specified or both 'SIGMA' and 'w' must be. 
%     There are a few exceptions to this requirement:
%       - When kappa = 0, neither of the variables 'K', 'SIGMA' and 'w' need to be specified.
%       - When kappa > 0, 'K' is left empty and 'W_bold' is specified, then only 'w' need to be
%         given as input ('SIGMA' need not be given).
%     
%   o If the 'T' parameter is left empty the 'tX' parameter must be set.
%
% The additional parameters are now explained one-by-one.
%
%   SIGMA:      [(kappa+1) x (kappa+1) real] The covariance matrix of the i.i.d. Gaussian vectors
%
%                       W_bold_i := (
%                                     W((i+1)/n) - W(i/n), 
%                                     int_{i/n}^{(i+1)/n} K((i+1)/n-s)dW(s),
%                                     int_{i/n}^{(i+1)/n} K((i+2)/n-s)dW(s),
%                                     ...,                                                  (4)
%                                     int_{i/n}^{(i+1)/n} K((i+kappa)/n-s)dW(s)
%                                   )'
%
%               for i = 0,1,...,floor(n*T)-1.
%
%               Remarks:
%                   o The parameter can be left empty (or unspecified) if instead the 'K'  
%                     parameter is set (or if the 'W_bold' parameter is used).
%                   o If kappa = 0 the vector in (4) should be interpreted as the 1x1 vector 
%                     containing just W((i+1)/n) - W(i/n).
%                   o If the SIGMA-matrix (either way inputted or computed) is not positive 
%                     definite (possibly due to round off errors) we use the modified Cholesky 
%                     factorisation algorithm from (Cheng and Higham, 1998) with the code obtained 
%                     from https://github.com/higham/modified-cholesky (retrieved on the 17th of 
%                     May 2020).
%
%   w:          [kappa x 1 real] A vector containing the integrals 
%
%                   w_k = int_{(k-1)/n}^{k/n} K(s)ds for k = 1,2,...,kappa
%
%               here assuming kappa > 0. If kappa = 0 it can be left empty. Can also be left  
%               empty (or unspecified) in general in which case the integrals are computed using 
%               numerical integration. In this case the 'K' parameter must be specified.
%
%   K:          [1x1 KernelFunctionClass] Kernel function. Must be vectorized. Only used if the 
%               'SIGMA' or 'w' parameter are left unspecified (or empty) in which case the 
%               SIGMA-matrix and/or w-vector will be computed by numerical integration. 
%
%   returnU:    [1x1 logical] This parameter should be set to true if the output should also 
%               contain simulated values of the U-factors defined as
%       
%                            U^j(t) =   int_0^t exp(-gamm_j*(t-s))b(s)ds 
%                                     + int_0^t exp(-gamm_j*(t-s))sigma(s)dW(s)
%
%               for j = 1,...,m.
%
%               The simulated values will be stored in the output variable 'Upaths'.
%
%               The default value for this parameter is false in which case the 'Upaths' output is 
%               left empty. If you do not need the U-factors it is recommended to keep 
%               returnU = false for better performance.
%
%   returndW:   [1x1 logical] If true, the output variable 'dW' is returned non-empty. If false 
%               (default) it is left empty. See also the description under 'outputs'. It is 
%               recommended (for performance) to leave returndW = false unless the output is 
%               needed. 
%
%   tX:         [1 x MX real] Parameter can only, and is exactly required to be used, if the 'T' 
%               parameter is left empty. The 'tX' parameter sets which time points to return  
%               simulated values of X(t) for. Together with the 'tU' parameter (possibly empty,  
%               see also below) we overwrite T = max([tX,tU]) and then use the hybrid multifactor 
%               scheme to first simulate the X-process on the grid points stated in (3). To return 
%               X-values at the time-points given by the 'tX' parameter, we then truncate the 'tX' 
%               time-points down to the nearest grid point in (3) and return the simulated X-values 
%               from those time points.
%
%               The vector must be sorted in ascending order and values must be unique.
%
%               If this parameter is instead unspecified or left empty (in which case the 'T' 
%               parameter is not) we instead return the X-values at the grid-points from (3).
%
%               See also the output variable 'Xpaths' for more information.
%
%   tU:         [1 x MU real] Parameter can only be used if returnU = true and the 'T' parameter  
%               is left empty. The 'tU' parameter sets which time points to return simulated values 
%               of the U-factors for. Together with the 'tX' parameter we overwrite T = max([tX,tU]) 
%               and then use the hybrid multifactor scheme to first simulate the processes on the 
%               grid points given by (3). To return values at the time points given by the 'tU' 
%               parameter, we then truncate the 'tU' values down to the nearest grid point in (3) 
%               and return the simulated values of the U-factors from those time points.
%
%               The vector must be sorted in ascending order and values must be unique.
%
%               If this parameter is instead left unspecified or empty (and returnU = true) the 
%               U-factors are returned at the grid points from (3).
%
%               See also the output variable 'Upaths' for more information.
%
%   Z:          [(M-1)*N x (kappa+1) real] i.i.d. standard normals needed for the simulation. 
%               Can be left unspecified or empty in which case they are automatically simulated in 
%               the code. If this parameter is used the W_bold-parameter (see below) must be 
%               left empty/unspecified.
%
%   W_bold:     [N x M-1 x kappa+1] Matrix containing the vectors in (4) and that for all the 
%               M-1 relevant time points and all paths. Should be used with care to ensure paths  
%               are simulated correctly. Inspect the code to see how to generate this matrix 
%               correctly. Parameter cannot be used simultaneously with the 'Z' parameter.
%
%   precision:  [1x1 string] The precision to perform the computations in. Options are 'double' 
%               and 'single'. Default is 'double'. To avoid unnecessary conversions it is 
%               recommended to input any of the following matrices in the same format as the 
%               specified precision: 
%               o b and sigma if given as Nx(M-1) matrices.
%               o Z or W_bold if inputted.
%
%	positive:	[1x1 logical] If true we replace any negative simulated X-values by its positive
%				part anywhere used in the scheme. Default is false.
%
%   explicit:   [1x1 logical] If true, we update the U-factors to a time point t_{i+1} as
%
%                   U^j(t_{i+1}) = U^j(t_i) + (b(t_i) - gamm_j*U^j(t_i))*dt 
%                                           + sigma(t_i)*(W(t_{i+1}) - W(t_i)).
%
%               for j=1,...,m.
%
%               If false, we update them as
%
%                   U^j(t_{i+1}) = (1/(1+gamm_j*dt))( U^j(t_i) + b(t_i)*dt 
%                                                              + sigma(t_i)*(W(t_{i+1})-W(t_i))
%                                                    )
%
%               for j=1,...,m.
%
%               Default is false.
%
% Output:
%   Xpaths: [1x1 struct] Simulated values of X(t). The struct has the following fields:
%
%               o t:            [1 x MX real] Time points. If only the 'T' parameter is specified 
%                               and not the 'tX' parameter, these time points are exactly those 
%                               shown in (3). Otherwise they are exactly the 'tX' time points.
%
%               o t_truncated:  [1 x MX real] Time points from which the corresponding values of X 
%                               are obtained. I.e. are the time points from the 't' field truncated 
%                               down to the nearest grid point. Note that if the 'T' parameter is 
%                               used, then fields t and t_truncated are equal.
%
%               o values:       [N x MX real] Simulated values of X(t).
%
%   Upaths: [1x1 struct] Simulated values of the U-factors. The struct has the fields to be 
%           explained below, but note that if returnU = false they are all left empty. If indeed
%           returnU = true we get:
%
%               o t:            [1 x MU real] Time points. If only the 'T' parameter is specified 
%                               and not the 'tU' parameter, these time points are exactly those 
%                               shown in (3). Otherwise they are exactly the 'tU' time points.
%
%               o t_truncated:  [1 x MU real] Time points from which the corresponding values of  
%                               the U-factors are obtained. Are the time points from the 't' field 
%                               truncated down to the nearest grid point. Note that if the 'T' 
%                               parameter is used, then fields t and t_truncated are equal.
%
%               o values:       [N x m x MU real] Simulated values of 
%       
%                                   U^j(t) =   int_0^t exp(-gamm_j*(t-s))b(s)ds 
%                                            + int_0^t exp(-gamm_j*(t-s))sigma(s)dW(s)
%
%                               for i = j,...,m.
%
%   dW:     [N x M-1 real or empty] If returndW = true this output contains the increments of the 
%           underlying Brownian motion W(t) across the grid points specified in (3). Otherwise it 
%           is left empty. 
%
% -------------------------------------------------------------------------------------------------
%  References:
% -------------------------------------------------------------------------------------------------
%   o Roemer, S.E.: Hybrid multifactor scheme for stochastic Volterra equations with completely
%     monotone kernels (2022). Working paper available at ssrn.com/abstract=3706253.
%   o Cheng, S.H., and Higham, N.J., A modified Cholesky algorithm based on a symmetric indefinite 
%     factorization. SIAM Journal on Matrix Analysis and Applications, 1998, 19(4), 1097-1110.
%
% -------------------------------------------------------------------------------------------------

% Parse name-value pair inputs:
p = inputParser;
addParameter(p,'SIGMA',[]);
addParameter(p,'w',[]);
addParameter(p,'K',[]);
addParameter(p,'returnU',false);
addParameter(p,'returndW',false);
addParameter(p,'tX',[]);
addParameter(p,'tU',[]);
addParameter(p,'precision','double',@(x) any(validatestring(x,{'single','double'})));
addParameter(p,'Z',[]);
addParameter(p,'W_bold',[]);
addParameter(p,'positive',false);
addParameter(p,'explicit',false);
parse(p,varargin{:});
v2struct(p.Results);

% Validate time point related inputs:
if ~isempty(T) && ~isempty(tX)
    error(['HybridMultifactorScheme: You cannot use both the ''T'' and ''tX'' ',...
          'parameter at the same time.']);
elseif isempty(T) && isempty(tX)
    error(['HybridMultifactorScheme: You must use either the ''T'' parameter or the ',...
          '''tX'' parameter.']);    
end

if ~isempty(tU) && ~isempty(T)
    error(['HybridMultifactorScheme: You cannot use both the ''T'' and ''tU'' ',...
          'parameter at the same time.']);    
end

if ~isempty(tU) && ~returnU
    error(['HybridMultifactorScheme: If the ''tU'' parameter is used',...
          ' the ''returnU'' parameter must be set to true.']);
end

if returnU && isempty(T) && isempty(tU)
    error(['HybridMultifactorScheme: If the returnU = true then either ',...
          'the ''T'' parameter or the ''tU'' parameter must be used.']);    
end

if ~isempty(tX)
    if any(diff(tX) <= 0)
        error(['HybridMultifactorScheme: The values of the ''tX'' vector must be strictly',...
               ' increasing.']);
    end
    if tX(1) < 0
        error('HybridMultifactorScheme: Negative time points are not allowed.');
    end
    if size(tX,1) > 1
        error('HybridMultifactorScheme: The ''tX'' parameter must be a row vector.');
    end
end
if ~isempty(tU)
    if any(diff(tU) <= 0)
        error(['HybridMultifactorScheme: The values of the ''tU'' vector must be strictly',...
               'increasing.']);
    end    
    if tU(1) < 0
        error('HybridMultifactorScheme: Negative time points are not allowed.');
    end    
    if size(tU,1) > 1
        error('HybridMultifactorScheme: The ''tU'' parameter must be a row vector.');
    end    
end

if isempty(T)
    T = max([tX,tU]);
end

% General set-up:
m = size(gamm,1);
dt = 1/n;
gamm = gamm';
c = c';

% The below computation avoids the round off errors that can arise using floor(n*T):
floor_nT = floor2(n*T);

if floor_nT <= 0
    error('HybridMultifactorScheme: Parameters must be such that floor(n*T) > 0.');
end

M = floor_nT + 1;
t_grid = (0:1:floor_nT)/n;

if size(tX,2) == 0
    tX = t_grid;
end
if size(tU,2) == 0 && returnU
    tU = t_grid;
elseif size(tU,2) == 0 && ~returnU
    tU = [];
end

MX = size(tX,2);
MU = size(tU,2);

% Adjust (if necessary) the requested time points to fit into the simulation grid. 
idxX = LocateTimepointsInDiscretisation(t_grid,tX);
if isempty(tU)
    idxU = [];
else
    idxU = LocateTimepointsInDiscretisation(t_grid,tU);
end
tX_trunc = t_grid(idxX);
tU_trunc = t_grid(idxU);

% Store the time points:
Xpaths.t = tX;
Upaths.t = tU;
Xpaths.t_truncated = tX_trunc;
Upaths.t_truncated = tU_trunc;

% Determine the input format of b:
if size(b,1) == N && size(b,2) == M - 1
    b_is_mat = true;
    b_is_function = false;
    b_is_constant = false;
elseif isa(b,'function_handle')
    b_is_mat = false;
    b_is_function = true;
    b_is_constant = false;
elseif size(b,1) == 1 && size(b,2) == 1
    b_is_mat = false;
    b_is_function = false;
    b_is_constant = true;
else
    error('HybridMultifactorScheme: ''b'' parameter is not valid.');
end

% Determine the input format of sigma:
if size(sigma,1) == N && size(sigma,2) == M - 1
    sigma_is_mat = true;
    sigma_is_function = false;
    sigma_is_constant = false;
elseif isa(sigma,'function_handle')
    sigma_is_mat = false;
    sigma_is_function = true;
    sigma_is_constant = false;
elseif size(sigma,1) == 1 && size(sigma,2) == 1
   sigma_is_mat = false;
   sigma_is_function = false;
   sigma_is_constant = true;
else
    error('HybridMultifactorScheme: ''sigma'' parameter is not valid.');
end

% Reformat g as a 1xM vector:
if size(g,1) == N && size(g,2) == M
    % do nothing
elseif isa(g,'function_handle')
    % overwrite g as a vector:
    g = g(t_grid);
elseif size(g,1) == 1 && size(g,2) == 1
    g = repmat(g,1,M);
else
    error('HybridMultifactorScheme: ''g'' parameter is not valid.');
end

if positive
	f_adj = @(x)(max(x,0));
else
	f_adj = @(x)(x);
end

if ~isempty(Z) && ~isempty(W_bold)
    error(['HybridMultifactorScheme: You cannot specify both the ''Z'' and ''W_bold'' ',...
           'parameter at the same time.']);
end

% Check kappa value is ok:
if mod(kappa,1) ~= 0
    error('HybridMultifactorScheme: Kappa must be integer valued.');
end
if kappa < 0 || kappa > floor_nT
    error('HybridMultifactorScheme: Kappa must be an integer between 0 and floor(n*T).');
end

% Initialize X-process:
Xpaths.values = zeros(N,MX,precision);
if idxX(1) == 1
    Xpaths.values(:,1) = g(:,1);
end

% Here we store temporary values of the X process:				   
X = g(:,1).*ones(N,1,precision);

% Initialize U-process:
if returnU
    Upaths.values = zeros(N,m,MU,precision);
else
    Upaths.values = [];
end

% Here we store temporary values of the U-factors:
U = zeros(N,m,precision);

% Sample random variables:
if isempty(W_bold)
    % Get the covariance matrix SIGMA:
    if isempty(SIGMA)
        if isempty(K) && kappa > 0
            error(['HybridMultifactorScheme: When the covariance matrix SIGMA is not ',...
                   'inputted, the W_bold parameter is unused, and kappa > 0, the kernel',...
                   'function K must be inputted instead.']);
        end
        SIGMA = ConvertMatrix(GetVolterraCovarianceMatrix({K},kappa,1/n),precision);

    else
        if size(SIGMA,1) ~= kappa + 1 || size(SIGMA,2) ~= kappa + 1
            error(['HybridMultifactorScheme: The covariance matrix SIGMA must be of the ',...
                   'size [(kappa + 1) x (kappa + 1)].']);
        end
        SIGMA = ConvertMatrix(SIGMA,precision);

    end    

    % Factorize covariance matrix: 
    % * An error will be thrown if SIGMA is not deemed positive semidefinite.
    A = chol_mod(SIGMA,10^(-14)); 

    % Sample random variables:
    if isempty(Z)
        Z = randn((M-1)*N,kappa+1,precision);
    else
        if size(Z,1) ~= (M-1)*N || size(Z,2) ~= kappa + 1
            error(['HybridMultifactorScheme: The Z-matrix containing the i.i.d. standard ',...
                   'normal random variables must be of the size [(M-1)*N x kappa + 1].']);
        end
        Z = ConvertMatrix(Z,precision);
    end
    W_bold = reshape(Z*(A.'),N,M-1,kappa+1);

else
    if size(W_bold,1) ~= N || size(W_bold,2) ~= M - 1 || size(W_bold,3) ~= kappa + 1
        error(['HybridMultifactorScheme: The pre-simulated W_bold-matrix has the ',...
               'wrong dimensions.']);
    end
    W_bold = ConvertMatrix(W_bold,precision);

end

% Initialize b and sigma matrices (if relevant):
if b_is_function
    b_mat = zeros(N,M-1,precision);
elseif b_is_mat
    b_mat = ConvertMatrix(b,precision);
end
if sigma_is_function
    sigma_mat = zeros(N,M-1,precision);
elseif sigma_is_mat
    sigma_mat = ConvertMatrix(sigma,precision);
end

% Get the w-vector:
if kappa > 0
    if isempty(w)
        if isempty(K)
            error(['HybridMultifactorScheme: When the w-vector is not inputted ',...
                   '(and kappa > 0) the kernel function K must be inputted instead.']);
        end
        w = ConvertMatrix(K.Integrate((0:kappa-1)'*dt,(1:kappa)'*dt,1,dt,c,gamm),precision).';        
    else
        if size(w,1) ~= kappa || size(w,2) ~= 1
            error('HybridMultifactorScheme: The w-vector must be of the size [kappa x 1].');
        end
        w = ConvertMatrix(w',precision);
    end
end

% Define some fixed vectors:
dummy1 = ConvertMatrix(c.*exp(-gamm*kappa*dt),precision);
dummy2 = ConvertMatrix((1./(1+gamm*dt)),precision);

% Evolve the X(t) process (and U-factors) forward in time:
% Remarks:
% o The code branches depending on b and/or sigma being scalars to achieve better performance.
% o We use the logical variable 'X_stored_in_Xpaths' to avoid copying/storing values of the X 
%   process unnecessarily. The code thus branches depending on whether the current iteration's  
%   X-values should be stored in the 'Xpaths' struct or the 'X' variable (or not at all).
% o In the first i=1,...,M-1 iterations we update U at t_(i-kappa} = (i-kappa)/n (when i > kappa)
%   and X at t_{i} = i/n. Note there that t_{M-1} = floor(n*T)/n is the last time point in the 
%   discretisation and thus these iterations allow us to obtain 'X' at any points in the 
%   discretisation. However, note also that these iterations only compute the U-factors at time 
%   points up to t_{floor(n*T)-kappa}. Thus, if returnU = true, we add another kappa iterations 
%   to compute U at the last few time points.
% o Assignments to the Upaths and Xpaths structs involve a call to Matlab's 'repmat' function in  
%   order to handle the case where multiple requested time points (inputted via either the tX or 
%   tU parameters) are truncated down to the same grid point. 
X_stored_in_Xpaths = false;
for i=1:(M-1+returnU*kappa)    
    % Update drift and diffusion coefficients to the time point t_{i-1} = (i-1)/n:
    if i <= M-1
        if b_is_function
            if X_stored_in_Xpaths
                % X-values from the last iteration was stored in the Xpaths struct:
                b_mat(:,i) = b(t_grid(i),f_adj(Xpaths.values(:,find(idxXi,1))));
            else
                b_mat(:,i) = b(t_grid(i),f_adj(X));
            end
        end
        if sigma_is_function
            if X_stored_in_Xpaths
                sigma_mat(:,i) = sigma(t_grid(i),f_adj(Xpaths.values(:,find(idxXi,1))));
            else
                sigma_mat(:,i) = sigma(t_grid(i),f_adj(X));
            end
        end
    end
    
    % Update the U-factors to the time point t_{i-kappa} = (i-kappa)/n:
    if i > kappa
        if ~b_is_constant && ~sigma_is_constant
            if explicit
                U = U.*(1 - gamm.*dt) + b_mat(:,i-kappa).*dt ...
                    + sigma_mat(:,i-kappa).*W_bold(:,i-kappa,1);
            else
                U = dummy2.*(U + (b_mat(:,i-kappa).*dt + sigma_mat(:,i-kappa).*W_bold(:,i-kappa,1)));
            end
        elseif b_is_constant && ~sigma_is_constant
            if explicit
                U = U.*(1 - gamm.*dt) + b.*dt + sigma_mat(:,i-kappa).*W_bold(:,i-kappa,1);
            else
                U = dummy2.*(U + (b.*dt + sigma_mat(:,i-kappa).*W_bold(:,i-kappa,1)));
            end
        elseif ~b_is_constant && sigma_is_constant
            if explicit
                U = U.*(1 - gamm.*dt) + b_mat(:,i-kappa).*dt + sigma.*W_bold(:,i-kappa,1);
            else
                U = dummy2.*(U + (b_mat(:,i-kappa).*dt + sigma.*W_bold(:,i-kappa,1)));
            end
        elseif b_is_constant && sigma_is_constant
            if explicit
                U = U.*(1 - gamm.*dt) + b.*dt + sigma.*W_bold(:,i-kappa,1);
            else
                U = dummy2.*(U + (b.*dt + sigma.*W_bold(:,i-kappa,1)));
            end
        end
        
        if returnU
            % Check if the newly simulated U-values should be stored.
            idxUi = ismember(idxU,i+1-kappa);
            if any(idxUi)
                Upaths.values(:,:,idxUi) = repmat(U,1,1,sum(idxUi));
            end
        end
    
    end
    
    % Update the X-process to the time point t_i:
    if i <= M-1
        % Check if we should store the value of X at time point t_i = i/n:
        idxXi = ismember(idxX,i+1);
        if ~any(idxXi) && ~b_is_function && ~sigma_is_function
            % We do not compute and store the X value at all. Continue to the next iteration:
            continue;
        elseif ~any(idxXi) && (b_is_function || sigma_is_function)
            % We compute and then store the value of the X-process in the 'X' object.
            X_stored_in_Xpaths = false;
        elseif any(idxXi)
            % We compute and then store the value of the X-process in the 'Xpaths' object.
            X_stored_in_Xpaths = true;
        end

        % Compute and store the value of X at the time point t_i:
        if kappa == 0
            % Pure lifted model: We only use the U-factors:
            if X_stored_in_Xpaths
                Xpaths.values(:,idxXi) = repmat(f_adj(g(:,i+1) + sum(c.*U,2)),1,1,sum(idxXi));
            else
                X = f_adj(g(:,i+1) + sum(c.*U,2));
            end
        else
            % Compute sum involving sigma(t,x):
            sigma_term = 0;
            for k=1:min(kappa,i)
                if ~sigma_is_constant
                    sigma_term = sigma_term + sigma_mat(:,i-k+1).*W_bold(:,i-k+1,1+k);
                else
                    sigma_term = sigma_term + sigma.*W_bold(:,i-k+1,1+k);
                end
            end
            kIdx = (1:min(kappa,i));
            % Add all the parts together. 
            % In general:  X(t) =   g(t)
            %                     + weighted sum of U-factors 
            %                     + sum involving b(t,x) values
            %                     + sum involving sigma(t,x) values
            % The U-factors are however absent for the first few iterations, see below.
            if i <= kappa
                % For the first few iterations the U-factors are not included:
                if ~b_is_constant
                    if X_stored_in_Xpaths
                        Xpaths.values(:,idxXi) = repmat(f_adj(g(:,i+1) ...
                                                           + sum(b_mat(:,i-kIdx+1).*w(kIdx),2) ...
                                                           + sigma_term),1,1,sum(idxXi));
                    else
                        X = f_adj(g(:,i+1) + sum(b_mat(:,i-kIdx+1).*w(kIdx),2) + sigma_term);
                    end
                else
                    if X_stored_in_Xpaths
                        Xpaths.values(:,idxXi) = repmat(f_adj(g(:,i+1) + b.*sum(w(kIdx),2) ...
                                                                 + sigma_term),...
                                                                 1,1,sum(idxXi));
                    else
                        X = f_adj(g(:,i+1) + b.*sum(w(kIdx),2) + sigma_term);
                    end
                end
            else
                % For the remaining iterations all terms are included:
                weighted_Us = dummy1.*U; % Matrix to store a temporary result (for performance)
                if ~b_is_constant
                    if X_stored_in_Xpaths
                        Xpaths.values(:,idxXi) = repmat(f_adj(g(:,i+1) + sum(weighted_Us,2) ...
                                                             + sum(b_mat(:,i-kIdx+1).*w(kIdx),2)...
                                                             + sigma_term),1,1,sum(idxXi));
                    else
                        X = f_adj(g(:,i+1) + sum(weighted_Us,2) ...
                                                 + sum(b_mat(:,i-kIdx+1).*w(kIdx),2)...
                                                 + sigma_term);
                    end
                else
                    if X_stored_in_Xpaths
                        Xpaths.values(:,idxXi) = repmat(f_adj(g(:,i+1) + sum(weighted_Us,2) ...
                                                                 + b.*sum(w(kIdx),2) ...
                                                                 + sigma_term),1,1,sum(idxXi));
                    else
                        X = f_adj(g(:,i+1) + sum(weighted_Us,2) ...
                                                 + b.*sum(w(kIdx),2) + sigma_term);
                    end
                end
            end
        end
    end
end

if returndW
    dW = squeeze(W_bold(:,:,1));
else
    dW = [];
end

end

