function [iv,se,S,optPrice,seOfPrice,idxCall] = GetOptionPrices(s0,y,q,K,T,N,n,V,dW1,S,tS,f,meanCV)
% Description: Implements Monte Carlo pricing for the model
%
%       dS(t) = S(t)(r(t) - delta(t))dt + S(t)sqrt(V(t))dW_1(t)
%
% where r(t) and delta(t) are deterministic functions (representing respectively the risk free 
% interest rate and the dividend rate) and V(t) a stochastic process.
%
% Main parameters:
%   s0:     [1x1 real] Initial asset price.
%   y:      [1x1 real or 1x1 CurveClass] Yield curve defined as  y(T) = (1/T)*int_0^T r(t)dt.
%           If set to a scalar the curve is assumed flat at that value.
%   q:      [1x1 real or 1x1 CurveClass] Dividend yield curve defined as 
%           q(T)=(1/T)*int_0^T delta(t)dt. If set to a scalar the curve is assumed flat at that 
%           value.
%   K:      [Nx1 real] Strikes.
%   T:      [Nx1 real] Expirations.
%   N:      [1x1 integer] The number of paths to use in the Monte Carlo estimator.
%   n:      [1x1 integer] The number of steps per year to simulate.
%   V:      [1x1 struct] Variance paths. 
%   dW1:    [...] Brownian motion driving the asset price.
%   S:      [... or empty] Asset price values already simulated. Can be left empty. Values need to
%           be normalized!
%   tS:     [...] Time points for S paths.
%   
% Output:
%   iv: [Nx1 real] Black-Scholes implied volatilities.
%   se: [Nx1 real] Standard errors (in volatility terms).
%

if ~exist('f','var');f = @(x)(x);end

nOpt = size(K,1);
[uniqT,idxUniqT] = unique(T);
maxT = max(uniqT);
dt = 1/n;

% Find grid point just below each expiration and define the grid points for the simulation:
t_grid_temp = (0:1:floor(n*maxT)+1)/n;
idxTadj = sum(T >= t_grid_temp,2);
Tadj = t_grid_temp(idxTadj).';
t_grid = t_grid_temp(1:max(idxTadj));

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
    error('GetOptionPrices: Expiries must be positive.');
end

% Validate inputs:
if size(K,2) > 1 || size(T,2) > 1 || size(K,1) ~= size(T,1)
    error('GetOptionPrices: Input contracts are not specified correctly.');
end

% Compute scaled strikes:
F = s0.*exp((y_eval-q_eval).*T);
ZCB = exp(-y_eval.*T);
K_adj = K./F;

% We estimate the time values by using out-of-the-money options only:
idxCall = K >= s0;
[optPrice,seOfPrice] = deal(zeros(nOpt,1));

% Simulate/compute S(t):
if ~exist('S','var') || isempty(S)
    S_inputted = false;
    dlogS = sqrt(V).*dW1 - 0.5*V*dt;
    S = zeros(N,size(t_grid,2));
    S(:, 1) = 1;
    S(:, 2:end) = exp(cumsum(dlogS,2));
else
    S_inputted = true;
end

for i=1:size(uniqT,1)
    idxT = T == uniqT(i);
    if ~S_inputted
        idxt = t_grid == Tadj(idxUniqT(i));
    else
        idxt = tS == Tadj(idxUniqT(i));
    end
    Ksub = K_adj(idxT);
    Fsub = F(idxT);
    ZCBsub = ZCB(idxT);
    idxCallSub = idxCall(idxT);
    
    X = zeros(sum(idxT),N);
    if any(idxCallSub)
        X(idxCallSub,:) = ZCBsub(idxCallSub).*Fsub(idxCallSub)...
                          .*max(S(:,idxt).'-Ksub(idxCallSub),0); 
    end
    if any(~idxCallSub)
        X(~idxCallSub,:) = ZCBsub(~idxCallSub).*Fsub(~idxCallSub)...
                          .*max(Ksub(~idxCallSub)-S(:,idxt).',0);     
    end
    
    if exist('meanCV','var')
        % Do regression:
        mu_X = mean(X.', 1); % Basic payoffs
        Y_cv = f(S(:,idxt));
        mu_Y = mean(Y_cv, 1); % Sample mean of CV
        diff_Y = Y_cv - mu_Y;
        alpha = -sum((X.' - mu_X).*diff_Y,1) ./ sum(diff_Y.^2,1);
        if any(isnan(alpha))
           alpha(isnan(alpha)) = 0;
        end 

        % Compute price estimates:
        X = ZCB(idxT).'.*F(idxT).'.*(X.' + (Y_cv - meanCV).*alpha);
        optPrice(idxT) = mean(X,1).'; 
        seOfPrice(idxT) = std(X,[],1).'./sqrt(size(X,1));
    else
        optPrice(idxT) = mean(X,2); 
        seOfPrice(idxT) = std(X,[],2)./sqrt(size(X,2));
    end
    
end

% Compute implied volatilities:
try
iv = blsimpv(s0,K,y_eval,T,optPrice,'Yield',q_eval,'Class',idxCall);


% Convert standard errors to implied volatilities:
se = seOfPrice./BSGreek('vega',[],s0,K,y_eval,T,iv,q_eval);
catch
    iv = NaN;
    se = NaN;
end
end

