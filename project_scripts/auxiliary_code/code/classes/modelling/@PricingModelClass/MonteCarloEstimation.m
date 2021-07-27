function [p,idxCall,se] = MonteCarloEstimation(~,k,ttm,cartProd,paths,optType,...
                                               condMC,cv,anti,prec,F,y,rho)
%
%   Estimates prices of call and put options given a set of (risk-neutral) paths.
%
%   Remark: The code assumes an asset price S(t) of the form
% 
%       dS(t) = S(t)(r(t) - delta(t))dt + S(t)sqrt(V(t))dW_2(t)               (*)
%
%   Here r(t) and delta(t) are deterministic functions (the risk-free interest rate and 
%   dividend yield respectively). The instantaneous variance V(t) must be adapted to the 
%   filtration generated some Brownian motion W_1(t) and the Brownian motion W_2(t) must 
%   then be given as 
%
%       W_2(t) = rho*W_1(t) + sqrt(1-rho^2)*W_perp(t)
%
%   where W_perp(t) is another Brownian motion that is independent of W_1(t) and -1 <= rho <= 1.
%
%   Price estimation is explained in detail below:
%  
%   The code uses many of the variance reduction methods from (McCrickerd and Pakkanen, 2018) 
%   where the rough Bergomi model is considered. The methods are however still valid in the 
%   more general setup explained above.
% 
%   Let K be the strike of an option and let w = 1 if it is a call option and w = -1 if it is 
%   a put option. Let us also write/define
%
%       y(t) = (1 / t) * int_0^t r(s) ds
%       q(t) = (1 / t) * int_0^t delta(s) ds.
%
%   Then under the stated assumptions we may rewrite the price as
%
%       exp(-y(T)*T)*E[(w*(S(T) - K))^+] = S(0)*exp(-q(T)*T)*E[(w*(tilde(S)(T) - exp(k)))^+]
%
%   with k = log(K/F(T)) being the log-moneyness and F(T) being the expiry T forward price 
%   and where tilde(S)(t) solves
%
%       d tilde(S)(t) = sqrt(V(t)) tilde(S)(t) dW_2(t), tilde(S)(0) = 1.
% 
%   It therefore suffices to perform the price estimation under the assumption of zero interest 
%   rate and dividend and an initial asset price of 1 and all that using the strike exp(k). 
%   The price can then easily be transformed to the more general setting by a simple 
%   multiplication. This is also how the code is structured.
%
%   In the following we therefore assume S(0) = 1 and r(t) = delta(t) = 0.
%
%   Prices are then estimated as
% 
%       price est. = average( X_i + alpha*(Y_i-E[Y])); i=1,...,# paths)      (**)
% 
%   where X is the main estimator of the option price, i.e. it will satisfy 
%   E[X] = E[{w*(S(T) - K)}^+]. See the parameter 'condMC' for the possible choices. Also, 
%   Y is a control variate (see the parameter 'cv' for the possible choices). The 'alpha' 
%   parameter is chosen in an asymptotically optimal way, see also (McCrickerd et al, 2018).
% 
% -----------------------------------------------------------------------------------------------    
% Parameters
% -----------------------------------------------------------------------------------------------    
% k:        [Nx1 real]       Log-moneyness.
%
% ttm:      [Mx1 real]       Expiries.
%
% cartProd: [1x1 logical]    If true then we return prices for the cartesian product of the k 
%                            and ttm vectors. Else (and here assuming N = M) we return prices 
%                            for each element (row-by-row).
%
% paths:     [1x1 struct]    Struct containing some of the following fields (We let A = # paths 
%                            including antithetic ones, B = # time points).
%
%                            o t:  [1xB real]    Time points. Must exactly include the unique  
%                                                values in the 'ttm' parameter and be sorted in 
%                                                ascending order.
%
%                            o S:  [AxB real]    The stochastic process S(t) as defined in the 
%                                                description.
%
%                            o S1: [AxB real]    The stochastic process S_1(t) defined as
%
%                                                    S_1(t) = exp{rho*[int_0^t sqrt(V(u)) dW_1(u)] 
%                                                                 - 0.5*rho^2*[int_0^t V(u) du]}.
%
%                            o QV: [AxB real]    The stochastic process QV(t) defined as
%
%                                                    QV(t) = int_0^t V(u) du.
%               
%                            o Y:  [AxB real]    Shared control variate. Exactly required if
%                                                cv = 'explicit'.
%
%                            o EY: [1xB real]    Known expectation of 'Y'. Exactly required if
%                                                cv = 'explicit'.
%
%                            Paths must be simulated under the assumption of an initial value of 
%                            1 and no drift, i.e. as in (*) except r(t)=delta(t)=0 and S(0)=1. 
%                            The code will then automatically correct for this. If there are 
%                            antithetic paths, the first half of the paths must be the original 
%                            paths and the second half the antithetic ones. The fields that
%                            we require depend on the choice of the other parameters. Specifically 
%                            we require the fields 'S1' and 'QV' if either condMC = true or 
%                            cv = 'timer_option'. We only require the field 'S' if condMC = false.
%
% optType:  [1x1 string]     The type of option to estimate prices for. Possible values are 'call' 
%                            (call options), 'put' (put options) and 'otm' (out-of-the-money 
%                            options). Moneyness is here measured in terms of forward-moneyness so 
%                            an option is at-the-money when strike = forward price.
%
% condMC:   [1x1 logical]    If false then X = {w*(S(T) - K)}^+ in (**), otherwise
%           
%                                X = E[(w*[S(T) - K])^+ | F_T^1 ] 
%                                  = BS((1-rho^2)*int_0^T V(u) du;S_1(T),k).
%
%                            where F_T^1 is the sigma algebra generated by the Brownian motion 
%                            W_1(t) on the time interval [0,T] and BS(x;y,z) is the Black-Scholes 
%                            price (put or call) with total variance x, current (spot) asset 
%                            price y and strike z.
%
% cv:       [1x1 string]     Control variate to be used. Options are
%
%                            o 'timer_option': Here we set
%
%                                     Y = BS(rho^2*[Q-QV(T)];S_1(T),k) 
%
%                              in (**) where 
%
%                                     Q = sup{ {QV(T)}_i , i=1,...,# paths}
%
%                              is the supremum of the quadratic variation of S(t) over all 
%                              simulated paths.
%
%                            o 'asset_price': Here we set Y = S_1(T) in (**) if condMC = true 
%                               and otherwise set Y = S(T).
%
%                            o 'none': Here we set Y = 0 in (**).
%
%                            o 'explicit': Here we use Y (and E[Y]) as given in the 'paths' input.
%
% anti:     [1x1 logical]    Set to true if half of the paths are antithetic, otherwise set it 
%                            to false. In the former case the first half of the paths in the 
%                            'paths' parameter must be the original paths and the second half 
%                            the antithetic ones. This is important for a correct computation 
%                            of standard errors.
%
% prec:     [1x1 string]     Precision. Options are 'double' and 'single'.
%
% F:        [Mx1 real]       Forward prices.
%
% y:        [1x1 CurveClass] Yield curve (interest rates).
%
% rho:      [1x1 real]       Correlation parameter. Is exactly required if condMC = true.
%
% -----------------------------------------------------------------------------------------------    
% Outputs
% -----------------------------------------------------------------------------------------------    
%   p:          [NxM or Nx1 real]      Option prices.
%   idxCall:    [NxM or Nx1 logical]   True ~ call, false ~ put. 
%   se:         [NxM or Nx1 real]      Standard errors.
%
% -----------------------------------------------------------------------------------------------    
% References
% -----------------------------------------------------------------------------------------------    
%   - McCrickerd, R. and Pakkanen, M.S., Turbocharging Monte Carlo pricing for the rough 
%     Bergomi model. Quantitative Finance, 2018, 18(11), 1877-1886.
%
% -----------------------------------------------------------------------------------------------    

    %% Validation and initialization:
    if condMC && isempty(rho)
        error(['PricingModelClass:MonteCarloEstimation: ''rho'' parameter is required ',...
               'when conditional Monte Carlo is used.']);
    end
    
    % Using cartProd = true with 1x1 dimensional inputs will collapse some of 
    % the matrix dimensions in the intermediate computations (which we do not
    % want). Thus this (valid) adjustment to the inputs:
    if size(k,1) == 1 && size(ttm,1) ~= 1
        k = repmat(k,size(ttm,1),1);
        cartProd = false;
    elseif size(ttm,1) == 1 && size(k,1) ~= 1
        ttm = repmat(ttm,size(k,1),1);
        cartProd = false;
    end          
    
    
    %% Compute adjusted payoffs:
    [Z,idxCall] = GetAdjustedPayoffs(k,ttm,cartProd,paths,optType,condMC,cv,anti,prec,rho);
    

    %% Compute prices and standard errors:
    if cartProd
        % Remark: Z is a matrix in 3 dimensions.
        yields = ConvertMatrix(y.Eval(ttm),prec);
        ZCBs = ConvertMatrix(exp(-yields.*ttm),prec);        

        N = size(Z,1);
        dummy = F.*ZCBs';        
        Zscaled = Z.*dummy;
        p = permute(squeeze(mean(Zscaled)),[2,1]);
        se = permute(squeeze(std(Zscaled)),[2,1])./sqrt(N);        

    else
        % Remark: Z is a cell array of matrices (one for each unique expiry)
        [p,se] = deal(zeros(size(k,1),1));
        uniqT = unique(ttm);
        for i=1:size(uniqT,1)
            idxT = ttm == uniqT(i);
            
            % Compute price (adjusted for dividend yield and interest rate):
            yields = ConvertMatrix(y.Eval(uniqT(i)),prec);
            ZCB = ConvertMatrix(exp(-yields.*uniqT(i)),prec);
        
            Zscaled = ZCB*F(i)*Z{i};
            N = size(Zscaled,1);
            p(idxT) = mean(Zscaled);
            se(idxT) = std(Zscaled)./sqrt(N);
        end
    end
    
end 


function [Z,idxCall,uniqT] = GetAdjustedPayoffs(k,ttm,cartProd,paths,optType,condMC,cv,anti,prec,...
                                                rho)
% 
%   An auxiliary function computing the adjusted payoffs from which price estimates (and standard 
%   errors can be computed).
% 

    %% Validation and initialization:
    ttPaths = paths.t;
    if ~issorted(ttPaths,'ascend')
        error(['PricingModelClass:MonteCarloEstimation: Time points in ''paths'' object ',...
               'must be sorted in ascending order.']);
    end
    uniqT = unique(ttm);
    if any(~ismember(uniqT,paths.t)) || any(~ismember(paths.t,uniqT))
        error(['PricingModelClass:MonteCarloEstimation: Time points in ''paths'' object ',...
               'must exactly be unique maturities.']);
    end
    
    if strcmpi(cv,'timer_option')
        Qn = max(paths.QV, [], 1);
    end
    
    if isfield(paths, 'S1')
        numPaths = size(paths.S1, 1);
    else
        numPaths = size(paths.S, 1);
    end
    
    if cartProd
        %% Cartesian product:
        k = ConvertMatrix(k,prec);
        k_grid = repmat(k, 1, size(ttm, 1));
        switch optType
            case 'call'
                idxCall = true(size(k_grid));
                idxCall_vec = true(size(k));
            case 'put'
                idxCall = false(size(k_grid));
                idxCall_vec = false(size(k));
            case 'otm'
                idxCall = k_grid>=0;
                idxCall_vec = k>=0;
        end

        numk = size(k,1);
        numT = size(ttm,1);
        if strcmpi(cv,'timer_option')
            totVar_Y = repmat((rho^2)*(Qn - paths.QV),1,1,1);
        end

        K = repmat(reshape(exp(k),1,1,numk),1,1);

        if strcmpi(cv,'timer_option')
            Y = nan(numPaths,numT,numk,prec);
            Y(:,:,idxCall_vec) =  bscall_3d(paths.S1,K(:,:,idxCall_vec),0,1,totVar_Y);
            Y(:,:,~idxCall_vec) = bsput_3d(paths.S1,K(:,:,~idxCall_vec),0,1,totVar_Y);
            clear totVar_Y;
        elseif strcmpi(cv,'asset_price')
            if condMC
                Y = paths.S1;
            else
                Y = paths.S;
            end
        elseif strcmpi(cv,'explicit')
                Y = paths.Y;
                EY = paths.EY;
        end
        
        X = zeros(numPaths,numT,numk,prec);
        
        if condMC
            totVar_X = repmat((1 - rho^2)*paths.QV,1,1,1);
            if any(idxCall_vec)
                 X(:,:,idxCall_vec) = bscall_3d(paths.S1,K(:,:,idxCall_vec),0,1,totVar_X);
            end
            if any(~idxCall_vec)
                 X(:,:,~idxCall_vec) = bsput_3d(paths.S1,K(:,:,~idxCall_vec),0,1,totVar_X);
            end
            clear totVar_X;
        else
            if any(idxCall_vec)
                 X(:,:,idxCall_vec) = max(paths.S - K(:,:,idxCall_vec),0);
            end
            if any(~idxCall_vec)
                 X(:,:,~idxCall_vec) = max(K(:,:,~idxCall_vec) - paths.S,0);
            end
        end

        if strcmpi(cv,'timer_option')
            EY = zeros(1,size(Qn,2),size(k,1),prec);
            EY(:,:,idxCall_vec) = bscall(1,squeeze(K(:,:,idxCall_vec))',0,1,Qn'*rho^2);
            EY(:,:,~idxCall_vec) = bsput(1,squeeze(K(:,:,~idxCall_vec))',0,1,Qn'*rho^2,0);
        elseif strcmpi(cv,'asset_price')
            EY = 1;
        end
        
        if anti
            % Average the original and antithetic paths:
            Nindep = size(X,1)/2;
            X = 0.5*(X(1:Nindep,:,:) + X(Nindep+1:end,:,:));
            if ~strcmpi(cv,'none')
                Y = 0.5*(Y(1:Nindep,:,:) + Y(Nindep+1:end,:,:));
            end
        end     
        
        mu_X = mean(X, 1);
        if ~strcmpi(cv,'none')
            mu_Y = mean(Y, 1);
            diff_Y = Y - mu_Y;
            alpha = -sum((X - mu_X).*diff_Y,1)./sum(diff_Y.^2,1);
            if any(any(isnan(alpha)))
                alpha(isnan(alpha)) = 0;
            end
            
            clear diff_Y;
            
        else
            alpha = 0;
            EY = 0;
            Y = 0;
        end
        
        % Compute adjusted payoffs:
        Z = X + alpha.*(Y - EY);
        
    elseif ~cartProd
        %% Non-cartesian product:
        k = ConvertMatrix(k,prec);
        K = exp(k);
        ttm = ConvertMatrix(ttm,prec);
        switch optType
            case 'call'
                idxCall = true(size(k));
            case 'put'
                idxCall = false(size(k));
            case 'otm'
                idxCall = k>=0;
        end
        uniqT = sort(unique(ttm));
        Z = cell(size(uniqT));
        for i=1:size(uniqT,1)
            idxT = ttm == uniqT(i);
            numStrikes = sum(idxT);
            if condMC || strcmpi(cv,'timer_option')
                s1 = paths.S1(:, i);
            end
            if ~condMC
                s = paths.S(:,i);
            end
                
            idxCallSub = idxCall(idxT);
            Ksub = K(idxT);

            if strcmpi(cv,'timer_option')
                totVar_Y = (rho^2)*( Qn(i) - paths.QV(:,i));
                Y = zeros(numPaths, numStrikes, prec);
                if any(idxCallSub)
                    Y(:,idxCallSub) = bscall(s1,Ksub(idxCallSub)',0,1,totVar_Y);
                end
                if any(~idxCallSub)
                    Y(:,~idxCallSub) = bsput(s1,Ksub(~idxCallSub)',0,1,totVar_Y,0);
                end
                clear totVar_Y;
            elseif strcmpi(cv,'asset_price')
                if condMC
                    Y = s1;
                else
                    Y = s;
                end
            elseif strcmpi(cv,'none')
                Y = 0;
                EY = 0;
                alpha = 0;
            elseif strcmpi(cv,'explicit')
                Y = paths.Y(:,i);
                EY = paths.EY(i);
            end
            
            X = zeros(numPaths, numStrikes, prec);
            if condMC == true
                totVar_X = (1 - rho^2)*paths.QV(:,i);
                if any(idxCallSub)
                    X(:,idxCallSub) = bscall(s1,Ksub(idxCallSub)',0,1,totVar_X);
                end
                if any(~idxCallSub)
                    X(:,~idxCallSub) = bsput(s1,Ksub(~idxCallSub)',0,1,totVar_X,0);
                end
            else
                if any(idxCallSub)
                    X(:,idxCallSub) = max(s-Ksub(idxCallSub)',0);
                end
                if any(~idxCallSub)
                    X(:,~idxCallSub) = max(Ksub(~idxCallSub)'-s,0);
                end
            end

            if strcmpi(cv,'timer_option')
                EY = zeros(1, numStrikes, prec);
                if any(idxCallSub)
                    EY(idxCallSub) = bscall(1, Ksub(idxCallSub)',0,1,Qn(i)*rho^2)';
                end
                if any(~idxCallSub)
                    EY(~idxCallSub) = bsput(1,Ksub(~idxCallSub)',0,1,Qn(i)*rho^2);
                end
            elseif strcmpi(cv,'asset_price')
                EY = 1;
            end
            
            if anti
                % Average the original and antithetic paths:
                Nindep = size(X,1)/2;
                X = 0.5*(X(1:Nindep,:) + X(Nindep+1:end,:));
                if ~strcmpi(cv,'none')
                    Y = 0.5*(Y(1:Nindep,:) + Y(Nindep+1:end,:));
                end
            end
            mu_X = mean(X, 1);
            if ~strcmpi(cv,'none')
                mu_Y = mean(Y, 1);
                diff_Y = Y - mu_Y;
                alpha = -sum((X - mu_X).*diff_Y,1) ./ sum(diff_Y.^2,1);
                if any(isnan(alpha))
                   alpha(isnan(alpha)) = 0;
                end
            end
            
            % Compute adjusted payoffs:
            Z{i} = X + bsxfun(@times,alpha,(Y - EY));
            
        end

    end


end















