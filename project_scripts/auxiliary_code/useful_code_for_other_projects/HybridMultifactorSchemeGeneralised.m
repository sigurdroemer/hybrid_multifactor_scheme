function [Xpaths,Upaths,dW,dZ] = HybridMultifactorSchemeGeneralised(N,n,T,g,b,sigma,rho,K,gamm,...
                                                                    c,kappa,varargin)
% 
%   Simulates the one-dimensional (but decoupled) stochastic Volterra equation
%
%      X(t) = g(t) + sum_{i=1}^d ( int_0^t K_i(t-s)b_i(s)ds + int_0^t K_i(t-s)sigma_i(s)dW_i(s) )
%
%   where g(t) is a deterministic function and we for i=1,...,d have that b_i(t) = b(t,X(t)) and 
%   sigma_i(t) = sigma_i(t,X(t)) are functions, W_i(t) a Brownian motion.
%
%   Defining W(t) := (W_1(t),...,W_d(t))', we assume cov(dW(t)) = rho*dt where rho is a dxd 
%   correlation matrix. 
%
%   To simulate the model we use a natural generalisation of the hybrid multifactor scheme of 
%   (Roemer, 2021). The K_i functions are therefore used exactly (i.e. without approximation) for 
%   a number of 'kappa' intervals near the origin and approximated by a sum of exponentials on the 
%   remainder of the domain. Specifically, letting 'dt' be the step size, we for t >= kappa*dt and 
%   all i=1,...,d, approximate
% 
%       K_i(t) approx(=) sum_{j=1}^{m_i} c_{i,j}*exp(-gamm_{i,j}*t)
%
%   where c_i = (c_{i,1},...,c_{i,m_i}) and gamm_i = (gamm_{i,1},...,gamm_{i,m_i}) are appropriate 
%   coefficients.
%
%   Remarks: 
%    o The correlated Brownian vector W(t) is constructed as W(t) := B*Z(t) where 
%      Z(t) := (Z_1(t),...,Z_d(t))' is a vector of independent Brownian motions and B a dxd matrix 
%      constructed such that B*transpose(B) = rho.
%
%    o Function has been optimised for speed. Code may thus be less readable in some places.
% 
% ------------------------------------------------------------------------------------------------
%   Main parameters
% ------------------------------------------------------------------------------------------------
%   N:     [1x1 integer]  Number of paths to generate.
%   n:     [1x1 integer]  Number of steps per unit of time, i.e. steps are of length dt := 1/n.
%   T:     [1x1 real ...  Final timepoint to simulate till. The process will then be simulated on 
%           or empty]     the M := floor(n*T) + 1 grid points stated below:
%                         
%                            0 < 0/n < 2/n < ... < floor(n*T)/n.                              (1)
% 
%                         The parameter can only be left empty if both the 'tX' and 'tU' parameters 
%                         are instead specified. See the description under 'Additional parameters' 
%                         for an explanation of how to use those parameters instead.
%   g:     [1x1 function] See main description.
%   b:     [1xd cell]     Each element represents b_i(t,x), i=1,...,d, from the main description  
%                         and can be one of the following types:
%
%                           o [1x1 function]:  A function of two variables representing b_i(t,x). 
%                                              Must be vectorised in the second argument.
%
%                           o [1x1 real]:      We interpret b_i(t,x) equal to this constant value. 
%                                              Recommended if you want b_i(t,x) constant as code 
%                                              has been further optimised for this case.
%
%                         Remark: While the solution of your particular SVE may be positive in 
%                         theory, this will not automatically hold for the numerical scheme. If 
%                         you specify an element of 'b' as a function you should therefore either 
%                         (1) check that the function returns a valid value even if the x input 
%                         is negative or (2) use the 'positive' parameter described further down.
%
%   sigma: [1xd cell]     As for the 'b' parameter but for the sigma_i(t,x), i=1,...,d, functions.
%                         See also the remark about positivity under 'b'.
%   rho:   [dxd real]     Correlation matrix for W-vector; see also the main description. Must be
%                         a positive semidefinite matrix with 1's along the diagonal. 
%   K:     [1xd cell]     Kernel functions. Each element is a [1x1 KernelFunctionClass]. 
%   gamm:  [1xd cell]     Exponents for the sum of exponentials approximations; See the main 
%                         description. Each element number i must be of the form [1 x m_i real].
%   c:     [1xd cell]     Coefficients for the sum of exponentials approximations; See the main 
%                         description. Each element number i must be of the form [1 x m_i real].
%   kappa: [1x1 integer]  Number of sub-integrals to approximate the K_i's exactly. Thus kappa = 0 
%                         corresponds to a pure sum of exponentials approximation of each 
%                         kernel function.
%
% ------------------------------------------------------------------------------------------------
%   Additional parameters
% ------------------------------------------------------------------------------------------------
%   Must be given in name-value pairs. Also, all parameters below are optional except for the
%   following special case:
%     
%     o If the 'T' parameter is left empty the 'tX' parameter must instead be set. Variables 
%       'T' and 'tX' cannot be used simultaneously.
%
%     o If the 'rho' parameter is left empty the 'B' parameter must instead be set. Variables 
%       'rho' and 'B' cannot be used simultaneously.
%
%   The additional parameters are now explained one-by-one.
%
%   B:          [dxL real]    Specifies the correlation structure for W(t); see the main 
%                             description. 
%
%   returnU:    [1x1 logical] This parameter should be set to true if the output should also 
%                             contain simulated values of the 'U-factors' which we define as
%       
%                               U_{i,j}(t) =   int_0^t exp(-gamm_{i,j}*(t-s))b_i(s)ds 
%
%                                            + int_0^t exp(-gamm_{i,j}*(t-s))sigma_i(s)dW_i(s)
%
%                             for i=1,...,d, j=1,...,m_i.
%
%                             Values will be stored in the output variable 'Upaths'.
%
%                             The default value for this parameter is false in which case the 
%                             'Upaths' output is left empty. If you do not need the U-factors 
%                             it is recommended to keep returnU = false for better performance.
%
%   returndW:   [1x1 logical] If true, the output variable 'dW' is returned non-empty containing
%                             all the increments of the W-process. If false (the default) it is 
%                             returned empty. See also the description under 'outputs'. It is 
%                             recommended (for performance) to leave returndW = false unless needed.
%
%   returndZ:   [1x1 logical] As parameter 'returndW' except for the Z-process.
%
%   tX:         [1 x MX real] Parameter can only, and is exactly required to be used, if the 'T' 
%                             parameter is left empty. The 'tX' parameter sets which timepoints 
%                             to return simulated values of X(t) for. Together with the 'tU' 
%                             parameter (possibly empty, see also below) we overwrite 
%                             T = max([tX,tU]) and then use the hybrid multifactor scheme to first 
%                             simulate the X-process on the grid points stated in (1). To return 
%                             X-values at the timepoints given by the 'tX' parameter, we then 
%                             truncate the 'tX' timepoints down to the nearest grid point in (1) 
%                             and return the simulated X-values from those timepoints.
%
%                             The vector must be sorted in ascending order and values must be 
%                             unique.
%
%                             If this parameter is instead unspecified or left empty (in which 
%                             case the 'T' parameter is not) we instead return the X-values at 
%                             the grid-points from (1).
%
%                             See also the output variable 'Xpaths' for more information.
%
%   tU:         [1 x MU real] Parameter can only be used if returnU = true and the 'T' parameter  
%                             is left empty. The 'tU' parameter sets which timepoints to return 
%                             simulated values of the U-factors for. Together with the 'tX' 
%                             parameter we overwrite T = max([tX,tU]) and then use the 
%                             hybrid multifactor scheme to first simulate the processes on the 
%                             grid points given by (1). To return values at the timepoints given 
%                             by the 'tU' parameter, we then truncate the 'tU' values down to the 
%                             nearest grid point in (1) and return the simulated values of the 
%                             U-factors from those timepoints.
%
%                             Vector must be sorted in ascending order and values must be unique.
%
%                             If this parameter is instead left unspecified or empty (and 
%                             returnU = true) the U-factors are returned at the grid points from 
%                             (1).
%
%                             See also the output variable 'Upaths' for more information.
%
%   Y:          [1xd cell]    Underlying i.i.d. Gaussians. Cannot be used simultaneously with the
%                             'W_bold' parameter. Can be left unspecified or empty in which case 
%                             numbers are automatically simulated in the code (unless of course 
%                             W_bold is used instead). 
%
%                             The i'th element of Y refers to the i.i.d. Gaussians needed to 
%                             simulate Z_i(t) and its related Volterra integrals. When 'rho' is 
%                             specified the i'th element must have the form [kappa*d+1 x N(M-1)...
%                             ...real].
%                             When 'B' is instead specified the i'th element must have the form
%                             [kappa*sum(B(:,i)~=0) + 1*any(B(:,i)~=0) x N(M-1) real]. Note that 
%                             fewer random numbers are (in general) needed in the second case as we 
%                             appropriately account for zero-entries of the B-matrix. We do not 
%                             account for this when 'rho' is used as computing the factorisation
%                             rho = B*transpose(B) numerically, will introduce non-zero values that
%                             in theory would be zero, in practise not. The required format of 
%                             Y would therefore be hard to predict for an end-user and for this
%                             reason we do not account for zero entries when 'rho' is used.
%
%                             In conclusion: If you for performance reasons want to avoid inputting 
%                             more numbers than needed, you should make sure to compute B outside 
%                             the function and ensure that the relevant entries are detectable 
%                             as zero using the check B(i,j) == 0. You should then run this 
%                             function with 'B' as the input.
%
%   W_bold:     [2xd cell]    Cell array containing the Brownian increments dW(t) as well as 
%                             needed Gaussian Volterra integrals of W(t) convolved with K.
%                             More precisely: For i=1,...,d, W_bold{1,i} must contain a 
%                             matrix of size [N x (M-1) real] representing W_i((j+1)/n) - W_i(j/n)
%                             for j=0,1,...,M-1. Also, for i=1,...,d, W_bold{2,i} must contain
%                             a matrix of size [N x (M-1) x kappa real] representing the following
%                             Gaussian integrals:
%
%                                     int_{j/n}^{(j+1)/n} K((j+1)/n-s)dW_i(s),
%                                     int_{j/n}^{(j+1)/n} K((j+2)/n-s)dW_i(s),
%                                     ...,                                                  
%                                     int_{j/n}^{(j+1)/n} K((j+kappa)/n-s)dW_i(s)
%
%                             for j = 0,1,...,M-1.
%
%                             Should be used with care to ensure paths are simulated with 
%                             the correct properties. Inspect the code to see how exactly this 
%                             input can be generated correctly.
%
%                             Parameter cannot be used simultaneously with the 'Y' parameter.
%
%   precision:  [1x1 string]  Precision to perform computations in. Options are 'double' and 
%                             'single'. Default is 'double'. To avoid unnecessary conversions it 
%                             is recommended to input 'Y' and 'W_bold' in the same format as the 
%                             requested precision: 
%
%	positive:	[1x1 logical] If true, we truncate any negative simulated X-values at zero. 
%                             Default is false.
%
%   explicit:   [1x1 logical] If true, we update the U-factors to the timepoint t_{k+1} as
%
%                               U_{i,j}(t_{k+1}) =    U_{i,j}(t_k) 
%                       
%                                                   + (b_i(t_k) - gamm_{i,j}*U_{i,j}(t_k))*dt 
%
%                                                   + sigma_i(t_k)*(W_i(t_{k+1}) - W_i(t_k)).
%
%                             for i=1,...,d, j=1,...,m_i.
%
%                             If false, we update them as
%
%                               U_{i,j}(t_{k+1}) = (1/(1+gamm_{i,j}*dt))( 
%                                                                U_{i,j}(t_k) + b_i(t_k)*dt 
%                                                              + sigma_i(t_k)*(W(t_{k+1})-W(t_k))
%                                                            )
%
%                             for i=1,...,d, j=1,...,m_i.
%
%                             The recommended choice is false (also the default) as the implicit 
%                             version of the scheme tends to be more numerically stable  
%
% ------------------------------------------------------------------------------------------------
%   Outputs
% ------------------------------------------------------------------------------------------------
%   Xpaths: [1x1 struct] Simulated values of X(t). The struct has the following fields:
%
%               o t:            [1 x MX real] Timepoints. If only the 'T' parameter is specified 
%                               and not the 'tX' parameter, these timepoints are exactly those 
%                               shown in (1). Otherwise they are exactly the 'tX' timepoints.
%
%               o t_truncated:  [1 x MX real] Timepoints from which the corresponding values of X 
%                               are obtained. I.e. are the timepoints from the 't' field truncated 
%                               down to the nearest grid point. Note that if the 'T' parameter is 
%                               used, then fields t and t_truncated are equal.
%
%               o values:       [N x MX real] Simulated values of X(t).
%
%   Upaths: [1xd cell or empty] Simulated values of U-factors. Is only non-empty if returnU = true.
%           Each element is then a [1x1 struct] representing the values of U_i, i=1,...,d. Struct 
%           number 'i' contains the following fields:
%
%               o t:            [1 x MX real] Timepoints. If only the 'T' parameter is specified 
%                               and not the 'tU' parameter, these timepoints are exactly those 
%                               shown in (1). Otherwise they are exactly the 'tU' timepoints.
%
%               o t_truncated:  [1 x MX real] Timepoints from which the corresponding values of  
%                               the U-factors are obtained. Are the timepoints from the 't' field 
%                               truncated down to the nearest grid point. Note that if the 'T' 
%                               parameter is used, then fields t and t_truncated are equal.
%
%               o values:       [N x m_i x MU real] Simulated values of 
%       
%                                   U_{i,j}(t) =   int_0^t exp(-gamm_{i,j}*(t-s))b_i(s)ds 
%
%                                                + int_0^t exp(-gamm_{i,j}*(t-s))sigma_i(s)dW_i(s)
%
%                               for j=1,...,m_i.
%
%   dW:     [1xd cell or empty] If returndW = true this ouput contains the increments of the 
%           Brownian motion W(t) across the grid points specified in (3). Each element will then 
%           be of the form [N x M-1 real].
%
% ------------------------------------------------------------------------------------------------
%   References
% ------------------------------------------------------------------------------------------------
%   o Roemer, S.E.: Hybrid multifactor scheme for stochastic Volterra equations, 2021,
%     Working paper available at ssrn.com/abstract=3706253.
%
%   o Cheng, S.H., and Higham, N.J., A modified Cholesky algorithm based on a symmetric indefinite 
%     factorization. SIAM Journal on Matrix Analysis and Applications, 1998, 19(4), 1097-1110.
%
% ------------------------------------------------------------------------------------------------

% Parse name-value pair inputs:
p = inputParser;
addParameter(p,'B',[]);
addParameter(p,'returnU',false);
addParameter(p,'returndW',false);
addParameter(p,'returndZ',false);
addParameter(p,'tX',[]);
addParameter(p,'tU',[]);
addParameter(p,'precision','double',@(x) any(validatestring(x,{'single','double'})));
addParameter(p,'Y',[]);
addParameter(p,'W_bold',[]);
addParameter(p,'positive',false);
addParameter(p,'explicit',false);
parse(p,varargin{:});
v2struct(p.Results);

% Check cell inputs are cells:
if ~iscell(b) || ~iscell(sigma) || ~iscell(K) || ~iscell(gamm) || ~iscell(c)
    error('HybridMultifactorScheme: Inputs b, sigma, rho, K, gamm, c must be cell arrays.');
end

% Validate timepoint related inputs:
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
        error('HybridMultifactorScheme: Negative timepoints are not allowed.');
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
        error('HybridMultifactorScheme: Negative timepoints are not allowed.');
    end    
    if size(tU,1) > 1
        error('HybridMultifactorScheme: The ''tU'' parameter must be a row vector.');
    end    
end

if isempty(T)
    T = max([tX,tU]);
end

% Validate inputs of dimension d:
d = size(b,2);
if size(b,1) > 1 || size(sigma,1) > 1 || size(K,1) > 1 || size(gamm,1) > 1 || size(c,1) > 1 ...
        || ~iscell(b) || ~iscell(sigma) || ~iscell(K) || ~iscell(gamm) || ~iscell(c) ...
        || size(b,2) ~= d || size(sigma,2) ~= d || size(K,2) ~= d || size(gamm,2) ~= d ...
        || size(c,2) ~= d
   error('HybridMultifactorScheme: Inputs b, sigma, K, gamm and c must be of the form [1xd cell].');
end

% Verify that elements of c and gamm are row vectors
% * We also compute the number of terms in each vector at the same time.
m = zeros(size(gamm));
for i=1:d
    m(i) = size(gamm{i},2);
    if size(gamm{i},1) > 1 || size(c{i},1) > 1 || size(c{i},2) ~= m(i)
        error(['HybridMultifactorScheme: Elements of ''c'' and ''gamm'' must (1) be row vectors ',...
               'and (2) of the same size.']);
    end
end

% Validate the correlation matrix:
if ~isempty(rho) && ~isempty(B)
    error('HybridMultifactorScheme: You cannot input both ''rho'' and ''B'' at the same time.');
    
elseif ~isempty(rho) % Here rho is inputted
    if size(rho,1) ~= d || size(rho,2) ~= d
        error('HybridMultifactorScheme: Correlation matrix ''rho'' must have dimensions dxd.');
    end
    if any(diag(rho) ~= 1)
        error('HybridMultifactorScheme: Correlation matrix ''rho'' must have 1''s on the diagonal.');
    end
    try
        B = chol_mod(rho); % Fails if matrix is not positive semidefinite.
    catch
        error(['HybridMultifactorScheme: Correlation matrix ''rho'' was not deemed ',...
               'positive semidefinite.']);
    end
    L = d;
    
else % Here B is inputted directly.
    
    % Check dimensions:
    if size(B,1) ~= d
        error('HybridMultifactorScheme: ''B'' must have ''d'' rows.');
    end
    
    % Check induced correlation matrix is valid:
    rho = B*B.';
    tol = 10^(-6);
    if any(abs(diag(rho)-1) > tol)
        error('HybridMultifactorScheme: B*transpose(B) must have 1''s on the diagonal.');
    end    
    
    try
        chol_mod(B*B.'); % Fails if B*B' is not positive semidefinite
    catch
        error('HybridMultifactorScheme: B*transpose(B) was not deemed positive semidefinite.');
    end    
    
    L = size(B,2);
    
end

% In the following we set up the simulation grid:
dt = 1/n;

% The below computation avoids the round off errors that can arise using floor(n*T):
floor_nT = floor(n*T+eps(n*T));

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

% Adjust (if necessary) the requested timepoints to fit into the simulation grid. 
idxX = LocateTimepointsInDiscretisation(t_grid,tX);
if isempty(tU)
    idxU = [];
else
    idxU = LocateTimepointsInDiscretisation(t_grid,tU);
end
tX_trunc = t_grid(idxX);
tU_trunc = t_grid(idxU);

% Store the timepoints:
Xpaths.t = tX;
Xpaths.t_truncated = tX_trunc;

[b_is_constant,sigma_is_constant] = deal(false(1,d));
for i=1:d
    % Determine the input format of the b's:
    if isa(b{i},'function_handle')
        b_is_constant(i) = false;
    elseif size(b{i},1) == 1 && size(b{i},2) == 1
        b_is_constant(i) = true;
    else
        error(['HybridMultifactorScheme: Format of ',num2str(i),...
               '''th ''b'' parameter is invalid.']);
    end

    % Determine the input format of the sigma's:
    if isa(sigma{i},'function_handle')
        sigma_is_constant(i) = false;
    elseif size(sigma{i},1) == 1 && size(sigma{i},2) == 1
        sigma_is_constant(i) = true;
    else
        error(['HybridMultifactorScheme: Format of ',num2str(i),...
               '''th ''sigma'' parameter is invalid.']);
    end

end

if positive
	f_adj = @(x)(max(x,0));
else
	f_adj = @(x)(x);
end

if ~isempty(Y) && ~isempty(W_bold)
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
    Xpaths.values(:,1) = g(0);
end

% Here we store temporary values of the X process:				   
X = g(0)*ones(N,1,precision);

% Initialize U-process:
if returnU 
    Upaths = cell(1,d);
    for i=1:d
        Upaths{i} = struct;
        Upaths{i}.t = tU;
        Upaths{i}.t_truncated = tU_trunc;        
        Upaths{i}.values = zeros(N,m(i),MU,precision);
    end
else
    Upaths = [];
end

% Here we store temporary values of the U-factors:
U = cell(1,d);
for i=1:d
    U{i} = zeros(N,m(i),precision);
end

% Sample Volterra integrals:
% * The purpose of the below code is to compute all needed Volterra sub-integrals of the form 
%   K_i*dW_i as well as Brownian increments dW_i (i=1,...,d). The results will be stored in
%   the cell array 'W_bold' and will then be ready to be used for the actual scheme.
%
% * The computations are separated into two steps. In the first step we compute for i=1,...,d 
%   and j=1,...,L all Volterra sub-integrals of the form K_i*dZ_j as well as Brownian increments 
%   dZ_j; the results of the first step will be stored in the object 'Z_bold'. In the second step
%   we aggregate the results of Z_bold using the B-matrix to obtain W_bold, i.e. to obtain the 
%   values of (K_i*dW_i,dW_i), i=1,...,d.
%
% * To more easily understand the code, recall also that all the mentioned Volterra integrals are 
%   needed at several shifted values of the kernel functions K (if kappa > 1). Note furthermore
%   that not all of the Volterra terms K_i*dZ_j may be needed. The code avoids the computation of 
%   such terms to improve performance.
if isempty(W_bold) 
    Z_bold = cell(d+1,L); % Each element will be a matrix containing values of dZ_j and K_i*dZ_j, 
                          % i=1,...,d, j=1,...,L. Values of the dZ_j's are stored in the first row. 
                          % Some elements representing the terms K_i*dZ_j may be left empty if the 
                          % corresponding B entry is empty.
                          
    W_bold = cell(2,d);   % The first row contains matrices for dW_i, i=1,...,d. The second row
                          % contains matrices for K_i*dW_i, i=1,...,d.
    
    % Step 1: We sample all sub-integrals of the form K_i*dZ_j (as well as dZ_j).
    for j=1:L 
        % * Samples are independent for different j (thus we loop over j)
        
        % First locate all K's we need to convolve with dZ_j:
        idxUse = find(B(:,j) ~= 0);

        if any(idxUse)
            % Construct covariance matrix:
            SIGMA = ConvertMatrix(GetVolterraCovarianceMatrix(K(idxUse),kappa,dt),precision);

            % Factorize covariance matrix:
            A = chol_mod(SIGMA);

            % Compute samples:
            if isempty(Y)
                dZ_and_KdZ = A*randn(size(SIGMA,1),N*(M-1),precision);
            else
                dZ_and_KdZ = A*ConvertMatrix(Y{j}(1:size(A,1),:),precision);
            end
            
            % Extract values:
            Z_bold{1,j} = reshape(dZ_and_KdZ(1,:,:),N,M-1);
            for i=1:size(idxUse,1)
                idxStart = 2 + (i-1)*kappa;
                idxEnd = 1 + i*kappa;
                Z_bold{idxUse(i)+1,j} = reshape(dZ_and_KdZ(idxStart:idxEnd,:),N,M-1,kappa); 
            end
            
        end
    end
    
    % Step 2: We compute values of K_i*dW_i and dW_i (for i=1,...,d):
    %  * Recall here that dW = B*dZ.
    for i=1:d
        for j=1:L
            if B(i,j) ~= 0
                if j==1
                    W_bold{1,i} = B(i,j)*Z_bold{1,j}; % ~ dZ_j
                    W_bold{2,i} = B(i,j)*Z_bold{i+1,j}; % ~ K_i*dZ_j                    
                else
                    W_bold{1,i} = W_bold{1,i} + B(i,j)*Z_bold{1,j}; % ~ dZ_j
                    W_bold{2,i} = W_bold{2,i} + B(i,j)*Z_bold{i+1,j}; % ~ K_i*dZ_j
                end
            end
        end
    end
    
end

% Initialize b and sigma matrices (if relevant):
b_mat = cell(1,d);
for i=1:d
    if ~b_is_constant(i)
        b_mat{i} = zeros(N,M-1,precision);
    end
end

sigma_mat = cell(1,d);
for i=1:d
    if ~sigma_is_constant(i)
        sigma_mat{i} = zeros(N,M-1,precision);
    end
end

% Get the w-vectors:
if kappa > 0
    w = cell(1,d);
    for i=1:d
        w{i} = ConvertMatrix(K{i}.Integrate((0:kappa-1)'*dt,(1:kappa)'*dt,...
                             1,dt,c{i},gamm{i}),precision);
    end
end

% Define some fixed vectors:
[dummy1,dummy2] = deal(cell(1,d));
for i=1:d
    dummy1{i} = ConvertMatrix(c{i}.*exp(-gamm{i}*kappa*dt),precision);
    dummy2{i} = ConvertMatrix((1./(1+gamm{i}*dt)),precision);
end

% Evolve the X-process (and U-factors) forward in time:
% Remarks:
% o The code branches depending on the b's and/or sigma's being scalars. We do this for better
%   better performance.
% o We use the logical variable 'X_stored_in_Xpaths' to avoid copying/storing values of the X 
%   process unnecessarily. The code thus branches depending on whether the current iteration's  
%   X-values should be stored in the 'Xpaths' struct or the 'X' variable (or not at all).
% o In the first i=1,...,M-1 iterations we update the U's at t_(i-kappa} = (i-kappa)/n (when i > 
%   kappa) and X at t_{i} = i/n. Note there that t_{M-1} = floor(n*T)/n is the last timepoint in the 
%   discretisation and thus these iterations allow us to obtain 'X' at any points in the 
%   discretisation. However, note also that these iterations only compute the U-factors at time 
%   points up to t_{floor(n*T)-kappa}. Thus, if returnU = true, we add another kappa iterations 
%   to compute U at the last few timepoints.
% o Assignments to the Upaths and Xpaths structs involve a call to Matlab's 'repmat' function in  
%   order to handle the case where multiple requested timepoints (inputted via either the tX or 
%   tU parameters) are truncated down to the same grid point. 
X_stored_in_Xpaths = false;
for i=1:(M-1+returnU*kappa)    
    % Update drift and diffusion coefficients to the timepoint t_{i-1} = (i-1)/n:
    if i <= M-1
        for j=1:d
            if ~b_is_constant(j)
                if X_stored_in_Xpaths
                    % X-values from the last iteration was stored in the Xpaths struct:
                    b_mat{j}(:,i) = b{j}(t_grid(i),f_adj(Xpaths.values(:,find(idxXi,1))));
                else
                    b_mat{j}(:,i) = b{j}(t_grid(i),f_adj(X));
                end
            end
            if ~sigma_is_constant(j)
                if X_stored_in_Xpaths
                    sigma_mat{j}(:,i) = sigma{j}(t_grid(i),f_adj(Xpaths.values(:,find(idxXi,1))));
                else
                    sigma_mat{j}(:,i) = sigma{j}(t_grid(i),f_adj(X));
                end
            end
        end
    end
    
    % Update the U-factors to the timepoint t_{i-kappa} = (i-kappa)/n:
    if i > kappa
        for j=1:d
            if ~b_is_constant(j) && ~sigma_is_constant(j)
                if explicit
                    U{j} = U{j}.*(1 - gamm{j}.*dt) + b_mat{j}(:,i-kappa).*dt ...
                        + sigma_mat{j}(:,i-kappa).*W_bold{1,j}(:,i-kappa);
                else
                    U{j} = dummy2{j}.*(U{j} + (b_mat{j}(:,i-kappa).*dt ...
                                         + sigma_mat{j}(:,i-kappa).*W_bold{1,j}(:,i-kappa)));
                end
            elseif b_is_constant(j) && ~sigma_is_constant(j)
                if explicit
                    U{j} = U{j}.*(1-gamm{j}.*dt) + b{j}.*dt ...
                                                 + sigma_mat{j}(:,i-kappa).*W_bold{1,j}(:,i-kappa);
                else
                    U{j} = dummy2{j}.*(U{j} + (b{j}.*dt ...
                                            + sigma_mat{j}(:,i-kappa).*W_bold{1,j}(:,i-kappa)));
                end
            elseif ~b_is_constant(j) && sigma_is_constant(j)
                if explicit
                    U{j} = U{j}.*(1 - gamm{j}.*dt) + b_mat{j}(:,i-kappa).*dt ...
                                                   + sigma.*W_bold{1,j}(:,i-kappa);
                else
                    U = dummy2{j}.*(U{j} + (b_mat{j}(:,i-kappa).*dt ...
                                         + sigma{j}.*W_bold{1,j}(:,i-kappa)));
                end
            elseif b_is_constant(j) && sigma_is_constant(j)
                if explicit
                    U{j} = U{j}.*(1 - gamm{j}.*dt) + b{j}.*dt + sigma{j}.*W_bold{1,j}(:,i-kappa);
                else
                    U{j} = dummy2{j}.*(U{j} + (b{j}.*dt + sigma{j}.*W_bold{1,j}(:,i-kappa)));
                end
            end

            if returnU
                % Check if the newly simulated U-values should be stored.
                idxUi = ismember(idxU,i+1-kappa);
                if any(idxUi)
                    Upaths{j}.values(:,:,idxUi) = repmat(U{j},1,1,sum(idxUi));
                end
            end
        
        end
    
    end
    
    % Update the X-process to the timepoint t_i:
    if i <= M-1
        % Check if we should store the value of X at timepoint t_i = i/n:
        idxXi = ismember(idxX,i+1);
        if ~any(idxXi) && all(b_is_constant) && all(sigma_is_constant)
            % We do not compute and store the X value at all. Continue to the next iteration:
            continue;
        elseif ~any(idxXi) && (any(~b_is_constant) || any(~sigma_is_constant))
            % We compute and then store the value of the X-process in the 'X' object.
            X_stored_in_Xpaths = false;
        elseif any(idxXi)
            % We compute and then store the value of the X-process in the 'Xpaths' object.
            X_stored_in_Xpaths = true;
        end

        % Compute and store the value of X at the timepoint t_i:
        if kappa == 0
            % Pure lifted model. We only use the U-factors:
            c_times_U_sum = 0;
            for j=1:d
                c_times_U_sum = c_times_U_sum + sum(c{j}.*U{j},2);
            end
            if X_stored_in_Xpaths
                Xpaths.values(:,idxXi) = repmat(f_adj(x0 + c_times_U_sum),1,1,sum(idxXi));
            else
                X = f_adj(x0 + c_times_U_sum);
            end
            
        else
            % Compute sum involving sigma(t,x):
            sigma_term = 0;
            for j=1:d
                for k=1:min(kappa,i)
                    if ~sigma_is_constant(j)
                        sigma_term = sigma_term + sigma_mat{j}(:,i-k+1).*W_bold{2,j}(:,i-k+1,k);
                    else
                        sigma_term = sigma_term + sigma{j}.*W_bold{2,j}(:,i-k+1,k);
                    end
                end
            end
            
            % Compute sum involving b(t,x):
            b_term = 0;
            for j=1:d
                for k=1:min(kappa,i)
                    if ~b_is_constant(j)
                        b_term = b_term + b_mat{j}(:,i-k+1).*w{j}(k);
                    else
                        b_term = b_term + b{j}.*w{j}(k);
                    end
                end
            end            
            
            % Compute weighted sum of U's:
            % * U's are though absent for the first few iterations (since kappa > 0).
            weighted_sum_of_Us = 0;
            if i > kappa
                for j=1:d
                    weighted_sum_of_Us = weighted_sum_of_Us + sum(dummy1{j}.*U{j},2);
                end
            end
                        
            % Add all the parts together. 
            % In general:  X(t) =   g(t)
            %                     + weighted sum of U-factors 
            %                     + sums involving b(t,x) values
            %                     + sums involving sigma(t,x) values
            if X_stored_in_Xpaths
                Xpaths.values(:,idxXi) = repmat(f_adj(g(t_grid(i)) ...
                                                      + weighted_sum_of_Us ...
                                                      + b_term ...
                                                      + sigma_term),...
                                                      1,1,sum(idxXi));
                                                  
            else
                X = f_adj(g(t_grid(i)) + weighted_sum_of_Us + b_term + sigma_term);
                
            end

        end
    end
end

if returndW
    dW = W_bold(1,1:end);
else
    dW = [];
end

if returndZ
    dZ = Z_bold(1,1:end);
else
    dZ = [];
end

end

