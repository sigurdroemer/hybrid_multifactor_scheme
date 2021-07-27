function [iv,se] = GetPricesRoughBergomi(s0,y,q,xi,eta,rho,H,K,T,scheme,N,n,kappa,varargin)
%
%   Implements the rough Bergomi model of (Bayer et al., 2016).
%
%   The model is briefly explained. Letting r(t) and delta(t) denote the deterministic risk-free 
%   interest rate and continuous proportional dividend yield (respectively) the asset price S(t) 
%   follows the dynamics
%
%       dS(t) = S(t)(r(t) - delta(t))dt + S(t)sqrt(V(t))dW_2(t)
%
%   under the risk-neutral measure. Here W_2(t) is a Brownian motion and V(t) the instantaneous 
%   variance process. The V(t) process is then modelled as
%
%       V(t) = xi(t)*exp(eta*sqrt(2H)*int_0^t(t-s)^{H-1/2}dW_1(s)-0.5*eta^2*t^{2H})
%
%   where W_1 is another Brownian motion s.t. dW_1(t) dW_2(t) = rho dt for a -1 <= rho <= 1. Also, 
%   xi(t) is a deterministic function, 0 < H < 1/2 and eta > 0.
%
% ------------------------------------------------------------------------------------------------
%  Main parameters
% ------------------------------------------------------------------------------------------------
%   s0:     [1x1 real] Initial asset price.
%
%   y:      [1x1 real or 1x1 CurveClass] Yield curve defined as  y(T) = (1/T)*int_0^T r(t) dt.
%           If set to a scalar the curve is assumed flat at that value.
%
%   q:      [1x1 real or 1x1 CurveClass] Dividend yield curve defined as 
%           q(T)=(1/T)*int_0^T delta(t) dt. If set to a scalar the curve is assumed flat at that 
%           value.
%
%   xi:     [1x1 real or 1x1 CurveClass] Forward variance curve. If set to a scalar value we assume
%           the curve flat at that value.
%
%   eta:    [1x1 real] Volatility-of-volatility parameter.
%
%   rho:    [1x1 real] Correlation parameter.
%
%   H:      [1x1 real] Hurst exponent.
%
%   K:      [Nx1 real] Strikes.
%
%   T:      [Nx1 real] Expirations.
%
%   scheme: [1x1 string]  Scheme specifying how to simulate the stochastic Volterra integral in 
%           the V(t) process. Options are 'HybridTBSS', we use the scheme from (Bennedsen et al., 
%           2017), and 'HybridMultifactor', we use the scheme from (Roemer, 2021).
%
%   N:      [1x1 integer] The number of paths to use in the Monte Carlo estimator. If turbo = true
%           this includes N/2 antithetic paths. In that case we require N to be divisible by 2.
%
%   n:      [1x1 integer] The number of steps per year to simulate.
%
%   kappa:  [1x1 integer] Corresponds to the 'kappa' from (Bennedsen et al., 2017) and (Roemer,
%           2021) depending on the scheme used.
%   
% ------------------------------------------------------------------------------------------------
%  Additional parameters (as name-value pairs)
% ------------------------------------------------------------------------------------------------
%   c:      [mx1 real] This parameter is only required if the 'HybridMultifactor' scheme is used.
%           Used as part of the exponential approximation where the kernel K(t) = t^(H-1/2) is
%           approximated by the function
%
%               K_hat(t) = sum_{i=1}^m c_i*exp(-gamm_i*t)                               (**)
%
%           where c = (c_1,...,c_m)' and gamm = (gamm_1,...,gamm_m)'.
%
%   gamm:   [mx1 real] This parameter is only required if the 'HybridMultifactor' scheme is used.
%           Used as part of the exponential approximation of the kernel function. See also (**).
%
%   Z1:     [floor(n*max(T))*N x kappa+1 or floor(n*max(T))*N/2 x kappa+1 real] i.i.d. standard 
%           normals needed for the simulation of the W_1(t) Brownian motion and thus the volatility 
%           process. Can be left unspecified in which case they are automatically simulated in the 
%           code. If turbo = false the input should be of the dimensions floor(n*max(T))*N x kappa+1
%           otherwise they should be floor(n*max(T))*N/2 x kappa+1.
%
%           Important remark: To avoid round-off errors one should compute the mathematical
%           expression floor(n*max(T)) as floor2(n*max(T)).
%
%   Z2:     [N x floor(n*max(T)) real] i.i.d. standard normals which generate a Brownian motion
%           W_perp(t) which is used to construct W_2(t) as:
%               W_2(t) = rho*W_1(t) + sqrt(1-rho^2)*W_perp(t). 
%           Can be left unspecified in which case the random variables are automatically simulated 
%           in the code. These random variables are only used/needed if turbo = false. 
%           Important remark: See the description of the 'Z1' parameter on how floor(n*max(T))
%           should be computed to avoid round-off errors.
%
%   turbo:  [1x1 logical] If true then we use the mixed estimator from (McCrickerd and Pakkanen, 
%           2018) to estimate prices, except we use 
%               
%                 S_1(t) = exp{rho*[int_0^t sqrt(V(u)) dW_1(u)] - 0.5*rho^2*[int_0^t V(u) du]}.
%
%           as the control variate instead of their timer-option based control variate. Note here 
%           that E(S_1(t)) = 1 for all t. Setting 'turbo' to true we also let half of the paths be
%           antithetic.
%
%           If the 'turbo' parameter is set to false, prices are estimated without any variance 
%           reduction techniques.
%
%           The parameter can be left unspecified in which case the code defaults to turbo = true.
%
%   SIGMA:  [(kappa+1) x (kappa + 1) real or empty] Covariance matrix to be used by the simulation 
%           scheme. Computed automatically if left empty (default), this is also the recommended
%           choice.
%
%   conv_method:
%           [1x1 logical or empty] Only applicable for scheme = 'HybridTBSS'. Specifies how the 
%           convolution should be computed. Options are 'fft' and 'conv_loop'. See the function
%           HybridTBSSScheme for more information. Uses 'fft' if left empty.
%
%   explicit:
%           [1x1 logical or empty] Scheme for the U-factors when using the hybrid multifactor 
%           scheme. Default is false. See the HybridMultifactorScheme function for more information.
%
% ------------------------------------------------------------------------------------------------
%  Output
% ------------------------------------------------------------------------------------------------
%   iv: [Nx1 real] Black-Scholes implied volatilities.
%   se: [Nx1 real] Standard errors (in volatility terms).
%
% ------------------------------------------------------------------------------------------------
%  References
% ------------------------------------------------------------------------------------------------
%   o Bayer, C., Friz, P. and Gatheral, J., Pricing under rough volatility. Quantitative Finance, 
%     2016, 16(6), 887-904.
%   o Bennedsen, M., Lunde, A. and Pakkanen, M.S., Hybrid scheme for Brownian semistationary 
%     procesess. Finance and Stochastics, 2017, 21(4), 931-965.
%   o McCrickerd, R. and Pakkanen, M.S., Turbocharging Monte Carlo pricing for 
%     the rough Bergomi model. Quantitative Finance, 2018, 18(11), 1877-1886.
%   o Roemer, S.E., Hybrid multifactor scheme for stochastic Volterra equations, 2021,
%     Working paper, available at https://www.ssrn.com/abstract=3706253.
%
% ------------------------------------------------------------------------------------------------

% Parse name-value pair inputs:
p = inputParser;
addParameter(p,'c',[]);
addParameter(p,'gamm',[]);
addParameter(p,'Z1',[]);
addParameter(p,'Z2',[]);
addParameter(p,'turbo',true);
addParameter(p,'SIGMA',[]);
addParameter(p,'conv_method',[]);
addParameter(p,'explicit',false);
parse(p,varargin{:});
v2struct(p.Results);

nOpt = size(K,1);
[uniqT,idxUniqT] = unique(T);
maxT = max(uniqT);
dt = 1/n;

% Find grid point just below each expiration and define the grid points for the simulation:
t_grid_temp = (0:1:floor2(n*maxT))/n;
idxTadj = sum(T >= t_grid_temp - eps(t_grid_temp),2);
Tadj = t_grid_temp(idxTadj).';
t_grid = t_grid_temp(1:max(idxTadj));
floor_nmaxT = max(idxTadj) - 1;

if isa(xi,'CurveClass')
    xi_eval = xi.Eval(t_grid(1:end-1));
else
    xi_eval = xi;
end
if isa(y,'CurveClass')
    y_eval = y.Eval(T);
else
    y_eval = y;
end
if isa(q,'CurveClass')
    q_eval = q.Eval(T);
else
    q_eval = q;
end

if min(T) < 0
    error('GetPricesRoughBergomi: Expiries must be positive.');
end

if turbo && mod(N,2) ~= 0
    error(['GetPricesRoughBergomi: If the ''turbo'' parameter is set to true the ''N'' ',...
           'parameter (the number of paths) must be divisible by 2.']);
end

% Validate inputs:
if xi < 0 || eta < 0 || abs(rho) > 1 || H <= 0.0 || H >= 0.5
    error('GetPricesRoughBergomi: One or more model parameters are invalid.'); 
end

if size(K,2) > 1 || size(T,2) > 1 || size(K,1) ~= size(T,1)
    error('GetPricesRoughBergomi: Input contracts are not specified correctly.');
end

% Compute scaled strikes:
F = s0.*exp((y_eval-q_eval).*T);
ZCB = exp(-y_eval.*T);
K_adj = K./F;

% We estimate the time values by using out-of-the-money options only:
idxCall = K >= s0;

if turbo
    N_indep = N/2;
else
    N_indep = N;
end

% Simulate variance process:
if strcmpi(scheme,'HybridTBSS')
    % Simulate the stochastic Volterra integral:
    [Y,~,dW1] = HybridTBSSScheme(N_indep,n,t_grid(end),sqrt(2*H)*eta,H-0.5,kappa,Z1,...
                                conv_method,[],[],SIGMA);
    
    % Compute the variance process:
    if ~turbo
        V = xi_eval.*exp(Y(:,1:end-1) - 0.5*t_grid(1:end-1).^(2.*H).*eta.^2);
    else
        V = zeros(N,floor_nmaxT);
        V(1:N/2,:) = xi_eval.*exp(Y(:,1:end-1) - 0.5*t_grid(1:end-1).^(2.*H).*eta.^2);
        V(N/2+1:end,:) = xi_eval.*exp(-Y(:,1:end-1) - 0.5*t_grid(1:end-1).^(2.*H).*eta.^2);
    end    
        
elseif strcmpi(scheme,'HybridMultifactor')
    % Compute covariance matrix and w-vector:
    if isempty(SIGMA)
        SIGMA = CovMatrixHybrid(n,kappa,H-0.5,'double');
    end
    if kappa > 0
        w = SIGMA(2:end,1);
    else
        w = [];
    end
    
    % Simulate the stochastic Volterra integral:
    [Y,~,dW1] = HybridMultifactorScheme(N_indep,n,t_grid(end),0,0,sqrt(2*H)*eta,gamm,c,kappa,...
                                        'Z',Z1,'returndW',true,'SIGMA',SIGMA,'w',w,'explicit',...
                                        explicit);
    
    % Compute the variance process:
    if ~turbo
        V = xi_eval.*exp(Y.values(:,1:end-1) - 0.5*t_grid(1:end-1).^(2.*H).*eta.^2);
    else
        V = zeros(N,floor_nmaxT);
        V(1:N/2,:) = xi_eval.*exp(Y.values(:,1:end-1) - 0.5*t_grid(1:end-1).^(2.*H).*eta.^2);
        V(N/2+1:end,:) = xi_eval.*exp(-Y.values(:,1:end-1) - 0.5*t_grid(1:end-1).^(2.*H).*eta.^2);
    end
    
else
    error('GetPricesRoughBergomi: The chosen scheme is not supported.');
end

[optPrice,seOfPrice] = deal(zeros(nOpt,1));
if turbo
    % We only need to simulate S_1(t):
    dlog_S1 = rho * sqrt(V) .* [dW1;-dW1] - 0.5 * (rho^2) * V * dt;
    S1 = zeros(N,size(t_grid,2));
    S1(:,1) = 1;
    S1(:,2:end) = exp(cumsum(dlog_S1, 2));    
    
    % Compute quadratic variation:
    QV = zeros(N,size(t_grid,2));
    QV(:,2:end) = cumsum(V,2).*dt;

    for i=1:size(uniqT,1)
        idxT = T == uniqT(i);
        idxt = t_grid == Tadj(idxUniqT(i));
        Ksub = K_adj(idxT);
        idxCallSub = idxCall(idxT);
        nOptSub = size(Ksub,1);
        
        % Compute conditional Monte Carlo estimator:
        totVar_X = (1 - rho^2)*QV(:,idxt);
        X_cond_mc = zeros(N,nOptSub);
        X_cond_mc(:,idxCallSub') = bscall(S1(:,idxt),Ksub(idxCallSub)',0,1,totVar_X);
        X_cond_mc(:,~idxCallSub') = bsput(S1(:,idxt),Ksub(~idxCallSub)',0,1,totVar_X);

        % Average the original and antithetic paths:
        X = 0.5*(X_cond_mc(1:N/2,:) + X_cond_mc(N/2+1:end,:));
        Y_cv = 0.5*(S1(1:N/2,idxt) + S1(N/2+1:end,idxt));
        
        % Do regression:
        mu_X = mean(X, 1);
        mu_Y = mean(Y_cv, 1);
        diff_Y = Y_cv - mu_Y;
        alpha = -sum((X - mu_X).*diff_Y,1) ./ sum(diff_Y.^2,1);
        if any(isnan(alpha))
           alpha(isnan(alpha)) = 0;
        end 

        % Compute price estimates:
        Z = ZCB(idxT).'.*F(idxT).'.*(X + (Y_cv - 1).*alpha);
        optPrice(idxT) = mean(Z).';
        
        % Compute standard errors:
        seOfPrice(idxT) = std(Z).'./sqrt(size(Z,1));
    end
else
    % Simulate/compute S(t):
    if isempty(Z1)
        dWperp = sqrt(dt).*randn(size(dW1));
    else
        if size(Z2,1) ~= N || size(Z2,2) ~= floor_nmaxT
           error('GetPricesrBergomi: The dimensions of the ''Z2'' parameter are invalid.'); 
        end
        dWperp = sqrt(dt).*Z2;
    end
    dW2 = rho.*dW1 + sqrt(1-rho^2).*dWperp;
    dlog_S = sqrt(V).*dW2 - 0.5*V*dt;
    S = zeros(N,size(t_grid,2));
    S(:, 1) = 1;
    S(:, 2:end) = exp(cumsum(dlog_S,2));
    
    for i=1:size(uniqT,1)
        idxT = T == uniqT(i);
        idxt = t_grid == Tadj(idxUniqT(i));
        Ksub = K_adj(idxT);
        Fsub = F(idxT);
        ZCBsub = ZCB(idxT);
        idxCallSub = idxCall(idxT);
        
        Z = zeros(sum(idxT),N);
        Z(idxCallSub,:) = ZCBsub(idxCallSub).*Fsub(idxCallSub)...
                          .*max(S(:,idxt).'-Ksub(idxCallSub),0); 
        Z(~idxCallSub,:) = ZCBsub(~idxCallSub).*Fsub(~idxCallSub)...
                          .*max(Ksub(~idxCallSub)-S(:,idxt).',0); 
        optPrice(idxT) = mean(Z,2); 
        seOfPrice(idxT) = std(Z,[],2)./sqrt(size(Z,2));
    end
end

% Compute implied volatilities:
iv = blsimpv(s0,K,y_eval,T,optPrice,'Yield',q_eval,'Class',idxCall);

% Convert standard errors to implied volatilities:
se = seOfPrice./BSGreek('vega',[],s0,K,y_eval,T,iv,q_eval);

end

