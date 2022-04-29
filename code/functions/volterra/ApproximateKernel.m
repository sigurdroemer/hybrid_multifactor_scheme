function [c,gamm,K,K_hat,exitflag] = ApproximateKernel(method,varargin)
%
% 	Approximates a completely monotone (kernel) function K by a sum of exponentials as
%
%       K_hat(t) = sum_j=1^m c_j*exp(-gamm_j*t)
% 
%   where (c_j, gamm_j), j=1,...,m, are coefficients.
%
%   To introduce notation let c = (c_1,...,c_m)' and gamm = (gamm_1,...,gamm_m)'.
%
%   Remarks: 
%    o Some methods are only implemented for the rough fractional kernel
%
%           K(t) = t^(H-0.5)/(Gamma(H + 0.5))                                           (*)
%
%      where H lies in (0,1/2).
%
%    o Some methods use numerical optimization. The function will only throw a warning and not an
%      error if something goes wrong in the optimization (e.g. the maximum number of iterations 
%      is reached). You should instead use the output 'exitflag' to check the optimization went 
%      well.
%
%    o It is generally recommended to use the method 'BM2005'.
%
% -------------------------------------------------------------------------------------------------
%  Parameters
% -------------------------------------------------------------------------------------------------
%  The parameter 'method' below sets the method to be used.
%
%   method: [1x1 string] The method for determining the (c,gamm) coefficients. Options:
%           
%               o 'AJEE2019': 
%                       We use the method of (Abi Jaber and El Euch, 2019) where the auxiliary mean
%                       reversion terms, denoted 'eta' in the paper, are set equidistantly. The
%                       spacing between the points is then chosen in a way that minimizes (an upper
%                       bound estimate) of the L2([0,T]) error. Method assumes that K is a rough
%                       fractional kernel.
%
%               o 'AJEE2019optim':
%                       We use the method of (Abi Jaber and El Euch, 2019) where the auxiliary mean
%                       reversion terms, denoted 'eta' in the paper, are numerically optimized 
%                       (without further restrictions) to minimize (an upper bound estimate) of 
%                       the L2([0,T]) error. Method assumes that K is a rough fractional kernel.
%
%               o 'AJ2019':
%                       We use the method of (Abi Jaber, 2019) where the r_m parameter (denoted r_n
%                       in that paper) is explicitly given as an input and the coefficients 
%                       (c, gamm) are then explicitly constructed from that. Method assumes that K 
%                       is a rough fractional kernel.
%
%               o 'AJ2019optim':
%                       We use the method of (Abi Jaber, 2019) where the r_m parameter (denoted r_n
%                       in that paper) is found by numerically minimizing the L2([delta,T]) error. 
%                       Method assumes that K is a rough fractional kernel.
%
%               o 'l2optim':
%                       We minimize the l2-error directly with respect to the coefficients (c,gamm)
%                       on a set of equidistant points over an interval of the form [delta,T].
%
%               o 'BM2005': 
%                       Uses the method of (Beylkin and Monzon, 2005) to locate a (nearly) minimal
%                       sum of exponentials approximation that approximately ensures an error below 
%                       some specified level. The error is computed across equidistant points in
%                       in an interval of the form [delta,T]. The error can either be measured
%                       as the normalised l2-error or as the relative uniform error.
%
% The remaining parameters must come in name-value pairs. Note that for a given method only some of
% them are required/used. We list them all below:
%
%   K:          [1x1 function] Kernel function. Must be completely monotone and vectorized. Can be
%               left unspecified if one wants to use the rough fractional kernel (*) and the 
%               parameter 'H' has instead been specified.
%
%   H:          [1x1 real] Fractional parameter. Must lie in the interval (0,1/2). Required for 
%               methods 'AJEE2019', 'AJEE2019optim', 'AJ2019', 'AJ2019optim'.
%
%   m:          [1x1 integer] Number of terms in the sum of exponentials approximation. 
%
%               Remark: Is not recommended for the method 'BM2005'. Use the 'epsilon' parameter
%               instead. For experimentation purposes the 'm' parameter can however be used even
%               for 'BM2005' (leave 'epsilon' unspecified in that case). Note though that the 
%               outputted number of terms may then be slightly smaller than the inputted 'm' value
%               if fewer than m roots are found in (0,1] for the m'th eigenpair.
%
%   error_measure:
%               [1x1 string] Only applies to the 'BM2005' method. Options are 'l2_normalised' and 
%               'uniform_relative'. With the choice 'l2_normalised' we seek a minimal representation 
%               that approximately ensures a normalised l2-error below 'epsilon'. With the choice 
%               'uniform_relative' we iterate until instead the relative uniform error is below 
%               'epsilon'. Parameter defaults to 'uniform_relative' if left unspecified.
%
%   epsilon:    [1x1 real] Only applies to the method 'BM2005'. Sets the error tolerance of the 
%               approximation. If error_measure = 'l2_normalised' we therefore seek a minimal 
%               representation ensuring ||K(t)-K_hat(t)||_2/||K(t)||_2 <= epsilon where ||.||_2 
%               denotes the euclidian norm and t is a vector of sampled points. If instead 
%               error_measure = 'uniform_relative' we seek a minimal representation ensuring 
%               max(abs(K(t)-K_hat(t))./K(t)) <= epsilon. Parameter cannot be used together with 
%               the 'm' parameter. If left empty and 'm' is left unspecified too then we default 
%               to epsilon = 10^(-3) if error_measure = 'l2_normalised' and epsilon = 10^(-2)
%               if error_measure = 'uniform_relative'.
%
%   n:          [1x1 integer] Only to be used by the methods 'l2optim' and 'BM2005'. Sets the number 
%               of intervals to divide [delta,T] into. A number n + 1 equidistant points will then 
%               be sampled. Remark: If the method is 'BM2005' and n is odd then we increment n by 
%               +1 to make it even.
%
%   T:          [1x1 real] The approximation is performed over the interval [delta,T] where 'delta' 
%               is a parameter explained below. This parameter is required for all methods.
%
%   delta:      [1x1 real] Kernel function will be approximated on the interval [delta,T].
%               Parameter is used by all other methods than 'AJEE2019', 'AJEE2019optim' and
%               'AJ2019'. Also, the methods 'l2optim' and 'BM2005' both require delta > 0.
%
%   rm:         [1x1 real] The r_m parameter from (Abi Jaber, 2019). Remark: m corresponds to n in
%               that paper. Only required for the method 'AJ2019'.
%
%   gamm0:      [mx1 real] Initial guess on the exponents gamm for the method 'l2optim'. Can be 
%               left empty in which case a default vector consisting of [0;0.1;0.2;...] is used.
%
%   nu:         [1x1 real] Only applies to the method 'BM2005' and only if the error_measure
%               parameter is set to 'uniform_relative'. Parameter specifies an initial guess on 
%               the normalised l2-error corresponding to a relative uniform error of epsilon. 
%               Defaults to epsilon/10 if left unspecified.
%
%   weight_tol: [1x1 real] Only applies to the method 'BM2005'. Let f denote K transformed to the 
%               domain [0,1] and let w_j, j=1,...,m, denote the weights of f (corresponding to
%               the c_j's for K). We then remove all weights where abs(w_j) < f(0)*weight_tol. 
%               Although, if all weights are deemed small per this inequality we do not remove the 
%               largest weight. Set to 0 to disable. Defaults to epsilon/10^4 if left unspecified.
%
%   approx_roots: 
%               [1x1 logical] Only used by the method 'BM2005'. If set to true then roots are found
%               by sampling Nz equidistant points in [0,1] and then running a local optimizer in  
%               each interval where we know a root exists. If set to false then we use Matlab's
%               'roots' function. Default is true if left unspecified.
%
%   Nz:         [1x1 integer] See the parameter 'approx_roots'. Default is 10^4 if left unspecified.
%
%   maxIter:    [1x1 logical] Only used by method 'BM2005' and only when the error measure is
%               set to 'uniform_relative'. Specifies the maximum number of iterations to find
%               a minimal representation ensuring a uniform relative error below epsilon.
%               Default is 20. 
%
%   options:
%               [1x1 struct] Options for fmincon optimizer when numerical optimization is used. 
%               Can be left empty in which case the default choice is:
%               options = optimset('Display','off');
%
%   disable_warn:   
%               [1x1 logical (optional)] For the methods that use numerical optimization we by
%               default throw a warning if something goes wrong during the optimization. 
%               See also the output parameter 'exitflag'. By setting this parameter to true, 
%               these warnings are disabled. Default value (if left unspecified) is false.
%
% -------------------------------------------------------------------------------------------------
%  Output
% -------------------------------------------------------------------------------------------------
%   c:        [mx1 real] Weights.
%   gamm:     [mx1 real] Exponents.
%   K:        [1x1 function] Kernel function.
%   K_hat:    [1x1 function] Approximating function.
%   exitflag: [1x1 integer] Exitflag from numerical optimization if such is used. 
%
% -------------------------------------------------------------------------------------------------
%  References
% -------------------------------------------------------------------------------------------------
%   o Abi Jaber, E., El Euch, O.: Multi-factor approximation of rough volatility models. SIAM 
%     Journal on Financial Mathematics, 2019, 10(2), 309-349.
%   o Abi Jaber, E.: Lifting the Heston model. Quantitative Finance, 2019, 19(12), 1995-2013.
%   o Beylkin, G., Monzon, L.: On approximation of functions by exponential sums. Applied and
%     Computational Harmonic Analysis, 2005, 19(1), 17-48.
%   o Roemer, S.E.: Hybrid multifactor scheme for stochastic Volterra equations (2022).
%     Working paper, available at ssrn.com/abstract=3706253.
%	   
% -------------------------------------------------------------------------------------------------
%

% Parse name-value pair inputs:
p = inputParser;
addParameter(p,'K',[]);
addParameter(p,'H',[]);
addParameter(p,'m',[]);
addParameter(p,'epsilon',[]);
addParameter(p,'n',[]);
addParameter(p,'T',[]);
addParameter(p,'delta',[]);
addParameter(p,'rm',[]);
addParameter(p,'options',[]);
addParameter(p,'disable_warn',false);
addParameter(p,'gamm0',[]);
addParameter(p,'Nz',[]);
addParameter(p,'approx_roots',[]);
addParameter(p,'weight_tol',[]);
addParameter(p,'error_measure',[]);
addParameter(p,'nu',[]);
addParameter(p,'maxIter',[]);
parse(p,varargin{:});
v2struct(p.Results);

if isempty(options)
    options = optimset('Display','off');
end

if isempty(K)
    if isempty(H)
        error('ApproximateKernel: If ''K'' is left unspecified then ''H'' should be set instead.');
    else
        K = @(t)(t.^(H - 0.5) ./ gamma(H + 0.5));
    end
end

exitflag = [];
switch method
    case 'AJEE2019'
        [c,gamm] = AbiJaberElEuch2019(H,T,m,false,options);
    case 'AJEE2019optim'
        [c,gamm,exitflag] = AbiJaberElEuch2019(H,T,m,true,options);
    case 'AJ2019'
        [c,gamm,exitflag] = AbiJaber2019(H,T,m,rm,delta,'explicit',options);
    case 'AJ2019optim'
        [c,gamm,exitflag] = AbiJaber2019(H,T,m,rm,delta,'optimize_L2',options);
    case 'l2optim'
        [c,gamm,exitflag] = Optimizel2(K,T,m,delta,n,options,gamm0);
    case 'BM2005'
        [c,gamm,exitflag] = BeylkinMonzon2005(K,delta,T,ceil(n/2),epsilon,m,approx_roots,Nz,...
                                              weight_tol,error_measure,nu,disable_warn,maxIter);
    otherwise
        error(['ApproximateKernel: Method ''', method, ''' is not supported.']);
end

if ~isempty(exitflag)
    if exitflag <= 0 && ~disable_warn
        warning(['ApproximateKernel: Optimizer returned exitflag ''', num2str(exitflag),'''. ',...
                 'It is recommended to check the solution.']);
    end
end

% Construct approximation as function:
if ~isempty(gamm) && ~isempty(c)
    [gamm,idxSort] = sort(gamm,'descend');
    c = c(idxSort);    
    K_hat = @(t)(sum(exp(-gamm'.*t).*c',2));
else
    K_hat = [];
end

end

function [c,gamm,exitflag,eta] = AbiJaberElEuch2019(H,T,m,optimL2,options)
% Description: Uses the methods of (Abi Jaber and El Euch, 2019) to approximate the rough fractional
% kernel.

exitflag = [];

% Validate inputs:
if H < 0 || H > 1/2
    error('ApproximateKernel:AbiJaberElEuch2019: ''H'' must be between 0 and 1/2.'); 
end
if T < 0
    error('ApproximateKernel:AbiJaberElEuch2019: ''T'' must be strictly positive.'); 
end
if m < 1
    error('ApproximateKernel:AbiJaberElEuch2019: ''m'' must be a strictly positive integer.'); 
end

% Define some useful functions:
dummy1 = gamma(H + 1/2)*gamma(3/2 - H);
dummy2 = gamma(H + 1/2)*gamma(1/2 - H)*(3/2 - H);
dummy3 = gamma(H + 1/2)*gamma(1/2 - H)*(5/2 - H);

mu_int = @(eta)((eta(2:end).^(0.5-H) - eta(1:end-1).^(0.5-H))./dummy1);
gamma_mu_int = @(eta)((eta(2:end).^(1.5-H) - eta(1:end-1).^(1.5-H)) ./ dummy2);
gamma_sqr_mu_int = @(eta)((eta(2:end).^(2.5-H) - eta(1:end-1).^(2.5-H)) ./ dummy3);

c_fun = @(eta)(mu_int(eta));
gamm_fun = @(eta)((1./c_fun(eta)).*gamma_mu_int(eta));

% Compute eta:
if ~optimL2
    pim = ((m^(-1/5))/T)*((sqrt(10)*(1-2*H))./(5-2*H))^(2/5);
    j = (0:m)';
    eta = j*pim;
    
else
    % Define some useful funtions:
    dummy4 = ((T.^(5/2))./(2*sqrt(5)));
    dummy5 = (H*gamma(H + 0.5)*gamma(0.5 - H)*sqrt(2)).^(-1);
    eta_from_delta = @(delta)([0;cumsum(delta)]);
    f2 = @(eta)(  dummy4 * sum( - 2*gamm_fun(eta).*gamma_mu_int(eta) ...
                             + (gamm_fun(eta).^2).*mu_int(eta) ...
                             + gamma_sqr_mu_int(eta) ...
                           ,1) ...
                   + dummy5*(eta(end))^(-H) );
            
    f2adj = @(delta)(f2(eta_from_delta(delta)));

    % Get initial guess and set bounds:
    [~,~,~,eta0] = AbiJaberElEuch2019(H,T,m,false,'off');
    delta0 = diff(eta0);
    lb = zeros(m,1);
    
    % Optimize and transform back to eta:
    [delta,~,exitflag] = fmincon(f2adj,delta0,[],[],[],[],lb,[],[],options);
    eta = eta_from_delta(delta);
    
end

% Transform to coefficients:
c = c_fun(eta);
gamm = gamm_fun(eta);

end

function [c,gamm,exitflag] = AbiJaber2019(H,T,m,rm,delta,type,options)
% Description: Uses the methods of (Abi Jaber, 2019) to approximate the rough fractional kernel.

exitflag = [];

% Validate inputs:
if H < 0 || H > 1/2
    error('ApproximateKernel:AbiJaberElEuch2019: ''H'' must be between 0 and 1/2.'); 
end
if T < 0
    error('ApproximateKernel:AbiJaberElEuch2019: ''T'' must be strictly positive.'); 
end
if m < 1
    error('ApproximateKernel:AbiJaberElEuch2019: ''m'' must be a strictly positive integer.'); 
end

if strcmpi(type,'explicit') && isempty(rm)
    error(['ApproximateKernel:AbiJaber2019: When type ''explicit'' is used the ''rm'' ',...
          'parameter must be set.']);
elseif strcmpi(type,'explicit') && rm <= 1
    error(['ApproximateKernel:AbiJaber2019: When type ''explicit'' is used the ''rm'' ',...
          'parameter must be strictly larger than 1.']);    
end

% Construct some useful functions:
alpha = H + 0.5;
j = (1:m)';
c_fun = @(rm)( (( (rm.^(1-alpha) - 1)*rm.^((alpha-1)*(1+m/2)) )./(gamma(alpha)*gamma(2-alpha)))...
            *  rm.^((1-alpha)*j) );
gamm_fun = @(rm)( ((1-alpha)/(2-alpha)) * ((rm^(2-alpha) - 1)/(rm^(1-alpha) - 1)) * rm.^(j-1-m/2)  );

K = @(t)(t.^(H - 0.5) ./ gamma(H + 0.5));
K_hat = @(t,rm)(sum(exp(-gamm_fun(rm).*t).*c_fun(rm),1));

% Find rm value:
if strcmpi(type,'optimize_L2')
    rm0 = 2.5;
    L2_err = @(rm)(integral(@(t)(( K(t) - K_hat(t,rm) ).^2),delta,T));
    [rm,~,exitflag] = fmincon(L2_err,rm0,[],[],[],[],1,[],[],options);
end

c = c_fun(rm);
gamm = gamm_fun(rm);

end

function [c,gamm,exitflag] = Optimizel2(K,T,m,delta,n,options,gamm0)
% Description: Numerically minimizes the l2-error to find the approximation.

% Define some useful functions and the initial guess:
K_hat = @(t,c,gamm)(sum(exp(-gamm.*t).*c,1));
dgamm_to_gamm = @(dgamm)(cumsum(dgamm));
if isempty(gamm0)
    dgamm0 = repmat(0.01,m,1);
else
    if size(gamm0,1) ~= m
        error('ApproximateKernel:Optimize: The initial guess ''gamm0'' must be of size mx1.');
    end
    gamm0 = sort(gamm0);
    dgamm0 = [gamm0(1);diff(gamm0)];
end

% Define error function:
dt = (T - delta)/n;
pts = (delta:dt:T);
K_pts = K(pts);
err_fun = @(par)(sum((K_pts-K_hat(pts,...
                      LeastSquaresWeightsFromGammas(K_pts,pts,dgamm_to_gamm(par)),...
                      dgamm_to_gamm(par))).^2));
lb = zeros(m,1);

% Optimize:
[par,~,exitflag] = fmincon(err_fun,dgamm0,[],[],[],[],lb,[],[],options);

% Extract solution:
gamm = dgamm_to_gamm(par);
c = LeastSquaresWeightsFromGammas(K_pts,pts,gamm);

end

function c = LeastSquaresWeightsFromGammas(K_pts,pts,gamm)
% Description: Used by the 'Optimize' function to find optimal weights when minimizing the l2-norm.
    A = exp(-pts'.*gamm');
    c = pinv(A)*K_pts.';
end

function [c,gamm,exitflag] = BeylkinMonzon2005(fOrig,delta,b,N,epsilon,m,approx_roots,Nz,...
                                               weight_tol,error_measure,mu,disable_warn,maxIter)
% Description: Implements the method of (Beylkin and Monzon, 2005).

% Initialize:
[c,gamm,exitflag] = deal([]);
if ~isempty(epsilon) && ~isempty(m)
    error(['ApproximateKernel:BeylkinMonzon2005: You cannot specify both ''epsilon'' and ''m''',...
          ' at the same time.']);
end
if isempty(N)
    error('ApproximateKernel:BeylkinMonzon2005: You must specify the ''N'' parameter.');
end
if isempty(error_measure);error_measure='uniform_relative';end
if isempty(epsilon)
    if strcmpi(error_measure,'l2_normalised')
        epsilon = 10^(-3);
    else
        epsilon = 10^(-2);
    end
end
if isempty(approx_roots);approx_roots=true;end
if isempty(Nz);Nz=10^4;end
if isempty(mu);mu=epsilon/10;end
if isempty(weight_tol);weight_tol=epsilon*10^(-4);end
if isempty(maxIter);maxIter=20;end

if ~isempty(m)
    force_one_iteration = true;
elseif strcmpi(error_measure,'l2_normalised')
    force_one_iteration = true;
    epsilon_l2 = epsilon;
elseif strcmpi(error_measure,'uniform_relative')
    force_one_iteration = false;
    epsilon_l2 = mu;
else
    error('ApproximateKernel:BeylkinMonzon2005: Unsupported error measure.');
end

if ~isempty(m)
    if m > N
        error('ApproximateKernel:BeylkinMonzon2005: We require m <= N ');
    end
    m_idx = m + 1;
end

% Grid for approximative root-finding:
if approx_roots
    z_grid = (0:(1/Nz):1)';
end

% Transform input function so it is defined on [0,1]:
f = @(x)(fOrig(x.*(b-delta)+delta));

% Compute h-vector:
ii = (0:2*N).';
x = ii./(2*N);
h = f(x);
h_norm = sqrt(sum(h.^2));

% Construct (N+1)x(N+1)-sized Hankel-matrix:
H = hankel(h(1:N+1),h(N+1:end));

% Compute only the first few eigenvalues and eigenvectors (for speed):
if isempty(m)
    num_eigs_init = 20;
else
    num_eigs_init = m + 1;
end
[V, Lambda, flag] = eigs(H,min(size(H,1),num_eigs_init));
if flag ~= 0
    error(['ApproximateKernel:BeylkinMonzon2005: Function call to ''eigs'' failed with ',...
             'flag = ', num2str(flag),'.']);
end
[sigmas,ind] = sort(max(diag(Lambda),0),'descend');

if isempty(m)
    % Guess the optimal m:
    m_idx = find(sigmas./h_norm <= epsilon_l2,1);
    m = m_idx - 1;
end

if isempty(m)
    % Compute the remaining eigenvalues and eigenvectors:
    [V, Lambda, flag] = eigs(H,size(H,1));
    if flag ~= 0
        error(['ApproximateKernel:BeylkinMonzon2005: Function call to ''eigs'' failed with ',...
               'flag = ', num2str(flag),'.']);
    end
    [sigmas,ind] = sort(max(diag(Lambda),0),'descend');

    % Guess again the optimal m:
    m_idx = find(sigmas./h_norm <= epsilon_l2,1);
    m = m_idx - 1;
    
    if isempty(m)
        m_idx = size(sigmas,1);
        m = m_idx - 1;
    end
end

% Iterate to find the minimal solution:
solution_found = false;
first_iter = true;
iter = 1;
while ~solution_found    
    % Find roots:
    u = flipud(V(:,ind(m_idx)));
    if ~approx_roots
        gamm_ = roots(u);

    else
        y = polyval(u,z_grid);
        pos = y > 0;
        sign_change = find(pos(1:end-1)~=pos(2:end));

        if isempty(sign_change)
            gamm_ = [];

        else        
            % Loop and locate roots:
            gamm_ = zeros(size(sign_change));
            fun = @(xIn)(polyval(u,xIn));
            for i=1:size(sign_change,1)
                [gamm_(i),~,exitflag] = fzero(fun,[z_grid(sign_change(i)),...
                                                   z_grid(sign_change(i)+1)]);
                if exitflag <= 0
                    error(['ApproximateKernel:BeylkinMonzon2005: Root finding failed with ',...
                            'exitflag = ', num2str(exitflag),'.']);
                end
            end

        end

    end
    gamm_ = unique(gamm_);

    % Filter roots:
    if ~isempty(gamm_)
        idxKeep = imag(gamm_) == 0 & real(gamm_) > 0 & real(gamm_) <= 1;
        gamm_ = gamm_(idxKeep);
        if ~isempty(gamm_)
            gamm_ = sort(gamm_,'descend');
            gamm_ = gamm_(1:min(m,size(gamm_,1)));
        end
    end

    if isempty(gamm_)
        gamm_ = 1;
    end

    % Determine weights by least-squares minimization:
    A = (gamm_.').^ii;
    w = pinv(A)*h;
    
    % Remove insignificant weights:
    idxRemove = abs(w) < f(0).*weight_tol;
    if any(idxRemove)
        if sum(idxRemove) == size(w,1)
            % Make sure at least one weight survives:
            [~,idxMaxRemove] = max(abs(w));
            idxRemove(idxMaxRemove) = false;
        end
        gamm_ = gamm_(~idxRemove);

        % Re-do least squares minimization:
        A = (gamm_.').^ii;
        w = pinv(A)*h;

    end    
        
    % Decide next step:
    if force_one_iteration
        % Stop with current solution:
        solution_found = true;
        
    else
        % Check the error:
        err = max(abs(h-A*w)./h);

        if (first_iter && err > epsilon) || (err > epsilon && err_prev > epsilon)
            % Here we attempt to increase m
            if m < N
                 m = m + 1;
                 m_idx = m_idx + 1;

                % Check if we need to increase the number of eigenvalues and eigenvectors:
                if m_idx > size(V,2)
                    % Find all eigenvalues and eigenvectors:
                    [V,Lambda,flag] = eigs(H,size(H,1));
                    [~,ind] = sort(max(diag(Lambda),0),'descend');
                    if flag ~= 0
                        error(['BeylkinMonzon2005: Function call to ''eigs'' failed with ',...
                                 'flag = ', num2str(flag),'.']);
                    end
                end

            else
                exitflag = -1;
                if ~disable_warn
                    warning(['ApproximateKernel: BeylkinMonzon2005: Algorithm was ',...
                             'unable to attain the desired precision.']);
                end
                return;

            end

        elseif (first_iter && err <= epsilon) || (err <= epsilon && err_prev <= epsilon)
            % Here we attempt to decrease m
            if m == 1
                solution_found = true;
            else
                m = m - 1;
                m_idx = m_idx - 1;
            end

        elseif err <= epsilon && err_prev > epsilon
            % Here we stop with the current solution
            solution_found = true;

        elseif err > epsilon && err_prev <= epsilon
            % Here we stop with the previous solution
            solution_found = true;
            gamm_ = gamm_prev;
            w = w_prev;

        end

        first_iter = false;
        
        % Save solution of current iteration:
        gamm_prev = gamm_;
        w_prev = w;
        err_prev = err;
        
    end
    
    iter = iter + 1;
    
    if iter > maxIter && ~solution_found
        if err <= epsilon
            solution_found = true;
            
        else
            exitflag = -1;
            if ~disable_warn
               warning(['ApproximateKernel: BeylkinMonzon2005: Algorithm was unable to attain ',...
                       'the desired precision within maxIter = ', num2str(maxIter),' iterations.']);
            end
            return;
            
        end
    end

end
        
% Transform solution to the domain [delta,b]:
t = 2*N*log(gamm_);
gamm = -t./(b-delta);
c = w.*exp(gamm.*delta);

end


