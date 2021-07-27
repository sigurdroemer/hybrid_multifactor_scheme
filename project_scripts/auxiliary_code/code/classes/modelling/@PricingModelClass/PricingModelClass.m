classdef PricingModelClass < GenericSuperClass
%
%    Superclass for pricing models. We assume an underlying asset with risk neutral dynamics
%    of the form
%
%       dS(t) = S(t)(r(t)-delta(t))dt + S(t)sqrt(V(t))dW(t)
%
%    where r(t) is the deterministic risk free interest rate, delta(t) the deterministic 
%    dividend yield, V(t) some appropriate stochastic process and W(t) a Brownian motion.
% 
% -----------------------------------------------------------------------------------------------
%  Properties
% -----------------------------------------------------------------------------------------------
%   s0:               [1x1 real]           Initial asset price.
% 
%   y:                [1x1 CurveClass]     Interest rate yield curve defined as
%                                               y(T) = (1/T) * int_0^T r(t) dt, T > 0.
%
%   q:                [1x1 CurveClass]     (integrated) dividend yield curve defined as
%                                               q(T) = (1/T) * int_0^T delta(t) dt, T > 0.
%
%   pricerSettings:   [1x1 Fourier...      Settings specifying which pricing algorithm should 
%                      PricerSettings...   be used and how. Which ones are available and how the 
%                      Class or ...        object should look will depend on the particular model. 
%                      1x1 MonteCarlo...
%                      PricerSettings...
%                      ... Class]     
%
    
properties
    s0
    y
    q
    pricerSettings
end
    
methods
function obj = PricingModelClass()
%
%    Constructor.
%
    if isempty(obj.y)
      obj.y = CurveClass('gridpoints',1,'values',0,'interpolation','linear','extrapolation','flat');
    end
    if isempty(obj.q)
      obj.q = CurveClass('gridpoints',1,'values',0,'interpolation','linear','extrapolation','flat');
    end
    if isempty(obj.s0)
        obj.s0 = 100;
    end

end
function [pS,pV,ivIdx] = GetPrices(obj,S_contracts,VIX_contracts,varargin)
%
%    Computes prices of puts and calls on S(t) and the VIX index. The latter is defined as
%
%       VIX(t) = sqrt((1/12) * int_t^{t+1/12} E[V(u)|F_t] du)
%
%    where F_t is the information set at time t.
%    
% -----------------------------------------------------------------------------------------------
%  Main parameters  
% -----------------------------------------------------------------------------------------------
%   * Remark: Exactly (and only) one of the fields (K, k, q) for parameters S_contracts and 
%     VIX_contracts (respectively) must be used.
%
%   S_contracts:       [1 x 1 struct]   Specification of S(t) options. The following fields apply:
%
%       o ttm:         [Ns x 1 real]    Maturities.
%       o K:           [Ms x 1 real]    Strikes.
%       o k:           [Ms x 1 real]    Log-moneyness, i.e. log(strike/forward).
%       o q:           [Ms x 1 real]    Strikes are chosen as quantiles of the terminal 
%                                       risk-neutral distribution of S(t). Values of 'q' must 
%                                       lie in the interval (0,1).
%
%       o cartProd:    [1x1 logical]    If 'true' we return prices for the cartesian product of 
%                                       the maturity and moneyness vectors. Otherwise we return
%                                       prices row-by-row (here assuming Ns = Ms). Field can be
%                                       left empty or unspecified in which case we default to 
%                                       'true' if Ns <> Ms and otherwise to 'false'.
%
%       o priceFormat: [1x1 string]     Format to return prices in. Options are 
%                                       'implied_volatility' (the default) and 'monetary'. The 
%                                       latter returns actual prices.
%
%       o optionType:  [1x1 string or   Type of options to return prices for. Options are 'call', 
%                       ... Ns x 1      'put', 'otm' (the default) and a [Ns x 1 logical] vector.
%                       ... logical]    The choice 'otm' returns prices for out-of-the-money 
%                                       options. If a logical vector is used we interpret
%                                       true ~ 'call', false ~ 'put'; in this case we also
%                                       require cartesian = false. Field is irrelevant if  
%                                       priceFormat = 'implied_volatility'.
%
%   VIX_contracts:     [1 x 1 struct]   Specification of VIX options. The following fields apply:
%       o ttm:         [Nv x 1 real]    Maturities.
%       o K:           [Mv x 1 real]    Strikes. 
%       o k:           [Mv x 1 real]    Log-moneyness, i.e. log(strike / vix futures).
%       o q:           [Mv x 1 real]    Strikes are chosen as quantiles of the terminal 
%                                       risk-neutral distribution of VIX(t).
%       o cartProd:    [1x1 logical]    Specifies how the maturity and moneyness inputs should
%                                       be interpreted. See 'S_contracts'.
%       o priceFormat: [1x1 string]     Price format, see 'S_contracts'.
%       o optionType:  [1x1 string or   Option type, see 'S_contracts'.
%                       ... Nv x 1
%                       ... logical]
%
%   Either of the parameters 'S_contracts' or 'VIX_contracts' can be left empty.
%
% -----------------------------------------------------------------------------------------------
%  Optional parameters  
% -----------------------------------------------------------------------------------------------
%  * Specify as name-value pairs
%
%   useExistingNumbers:     [1x1 logical]     Algorithm attempts to (re)use any existing random
%                                             numbers stored in the obj.pricerSettings.random
%                                             property instead of sampling new ones. Generally 
%                                             not intended for an end-user. Should be used
%                                             with care. Default is false. 
%
%   specialRun:             [1x1 logical]     Not intended for an end user. Used internally
%                                             by algorithm for repeated calls. Default is false.
%
%   dispIter:               [1x1 logical]     Shows a text in command window each time an iteration
%                                             is completed when obj.pricerSettings.numRepeat > 1.
%           
%
% -----------------------------------------------------------------------------------------------    
%  Outputs
% -----------------------------------------------------------------------------------------------    
%   pS: [1x1 struct]   Prices of options on S(t). Contains the following fields:
%      o p:            [Ns x 1 or Ns x Ms real]    Prices or implied volatility.
%      o priceFormat:  [1x1 string]                Format of price.
%      o idxCall:      [Ns x 1 or Ns x Ms logical] Option type. True ~ call, false ~ put.
%      o F:            [Ns x 1 or Ns x Ms real]    Forward prices.
%      o K:            [Ns x 1 or Ns x Ms real]    Strikes.
%      o k:            [Ns x 1 or Ns x Ms real]    Log-moneyness values.
%      o q:            [Ns x 1 or Ns x Ms real]    Requested quantiles (if relevant).
%      o ttm:          [Ns x 1 or Ns x Ms real]    Expirations.
%      o se:           [Ns x 1 or Ns x Ms real]    Standard errors.
%
%   pV: [1x1 struct] Prices of options on VIX(t). Contains the following fields:
%      o p:            [Nv x 1 or Nv x Mv real]    Prices or implied volatility.
%      o priceFormat:  [1x1 string]                Format of price.
%      o idxCall:      [Nv x 1 or Nv x Mv logical] Option type. True ~ call, false ~ put.
%      o F:            [Nv x 1 or Nv x Mv real]    Futures prices.
%      o K:            [Nv x 1 or Nv x Mv real]    Strikes.
%      o k:            [Nv x 1 or Nv x Mv real]    Log-moneyness values.
%      o q:            [Nv x 1 or Nv x Mv real]    Requested quantiles (if relevant).
%      o ttm:          [Nv x 1 or Nv x Mv real]    Expirations.
%      o se:           [Nv x 1 or Nv x Mv real]    Standard errors of option prices.
%      o se_F:         [Nv x 1 or Nv x Mv real]    Standard errors of futures prices.
%
%   ivIdx: [1x1 logical] Not intended for an end user. Used internally to handle repeated calls 
%                        to function.
% -----------------------------------------------------------------------------------------------    

   % Initialize outputs:
   [pS,pV,ivIdx] = deal(struct);
   
   % Determine which assets to price options on:
   if ~isempty(S_contracts);doPriceS=true;else;doPriceS=false;end
   if ~isempty(VIX_contracts);doPriceVIX=true;else;doPriceVIX=false;end
   if ~doPriceS && ~doPriceVIX;return;end

   % Parse variable inputs:
   p = inputParser;
   addParameter(p,'useExistingNumbers',false);
   addParameter(p,'specialRun',false);
   addParameter(p,'dispIter',false);
   parse(p,varargin{:});
   v2struct(p.Results);

   %% Parse and validate contract inputs:
   % Validate that no unexpected fields are specified:
   if ~isempty(S_contracts)
       validFields = {'ttm','k','K','q','cartProd','priceFormat','optionType'};
       actualFields = fieldnames(S_contracts);
       if any(~ismember(actualFields,validFields))
            error('PricingModelClass:GetPrices: Invalid field names in ''S_contracts'' input.');
       end
   end
   
   % Validate that no unexpected fields are specified:
   if ~isempty(VIX_contracts)
       validFields = {'ttm','K','k','q','cartProd','priceFormat','optionType'};
       actualFields = fieldnames(VIX_contracts);
       if any(~ismember(actualFields,validFields))
            error('PricingModelClass:GetPrices: Invalid field names in ''VIX_contracts'' input.');
       end
   end   
   
   % Basic check on dimensions of contract inputs:
   if (~isempty(S_contracts) && ~VerifyStructFieldDimensions(S_contracts,{'ttm','k','K','q'},...
       {'column_vector','column_vector','column_vector','column_vector'},...
       [false,true,true,true],[false,true,true,true])) ...
       ||(~isempty(VIX_contracts) && ~VerifyStructFieldDimensions(VIX_contracts,...
           {'ttm','k','K','q'},{'column_vector','column_vector','column_vector','column_vector'},...
           [false,true,true,true],[false,true,true,true]))
       
       error(['PricingModelClass:GetPrices: Contracts are incorrect specified. Please ',...
              'check that the specification of the contract inputs is valid.']);
          
   end
   
   % Additional checks:
   ttmUniq = [];
   if doPriceS
       if ~isfield(S_contracts,'priceFormat')
           S_contracts.priceFormat = 'implied_volatility';
       end
       q_is_used_S = isFieldNonEmpty(S_contracts,'q');
       k_is_used_S = isFieldNonEmpty(S_contracts,'k');
       K_is_used_S = isFieldNonEmpty(S_contracts,'K');
       if q_is_used_S + k_is_used_S + K_is_used_S ~= 1
            error(['PricingModelClass: Exactly one of the fields ''k'', ''K'' or ''q'' must be ',...
                   'specified for the ''S_contracts'' input.']);
       end
       if q_is_used_S
            moneyness_S = S_contracts.q;
            if any(S_contracts.q <= 0 | S_contracts.q >= 1)
               error('PricingModelClass:GetPrices: Quantiles must lie in (0,1).');
            end
       elseif k_is_used_S
            moneyness_S = S_contracts.k;
       elseif K_is_used_S
            moneyness_S = S_contracts.K;
       end
       if isFieldEmpty(S_contracts,'cartProd')
           if size(moneyness_S,1) ~= size(S_contracts.ttm,1)
               S_contracts.cartProd = true;
           else
               S_contracts.cartProd = false;
           end
       end
       if ~S_contracts.cartProd 
           if (k_is_used_S && ~(size(S_contracts.k,1) == size(S_contracts.ttm,1))) || ...
              (q_is_used_S && ~(size(S_contracts.q,1) == size(S_contracts.ttm,1))) || ...
              (K_is_used_S && ~(size(S_contracts.K,1) == size(S_contracts.ttm,1)))
               error(['PricingModelClass:GetPrices: When not using a Cartesian',...
                      ' product the moneyness and maturity vectors must be', ...
                      ' of the same length.']);   
           end
       end
       if S_contracts.cartProd && (any(size(moneyness_S) ~= size(unique(moneyness_S))) ||...
               any(size(S_contracts.ttm) ~= size(unique(S_contracts.ttm))))
                error(['PricingModelClass:GetPrices: The moneyness and expiry ',...
                       'vectors must be unique when a Cartesian product is ',...
                       'requested.']);     
       end
       if isFieldEmpty(S_contracts,'optionType')
           S_contracts.optionType = 'otm';
       elseif islogical(S_contracts.optionType)
           if S_contracts.cartProd
               error(['PricingModelClass:GetPrices: If the optionType ',...
                      'parameter is specified as a logical vector we must have ',...
                      'cartProd = false.']);
           elseif size(moneyness_S,1) ~= size(S_contracts.optionType,1) || ...
                   size(moneyness_S,1) ~= size(S_contracts.ttm,1)
               error(['PricingModelClass:GetPrices: If the optionType parameter is specified ',...
                      'as a logical vector, the optionType, moneyness and maturity vectors ',...
                      'must have the same sizes.']);
           end
       else
           if ~ismember(S_contracts.optionType,{'call','put','otm'})
               error(['PricingModelClass:GetPrices: If the parameter ',...
                      'optionType is specified as a string only the values ',...
                      '''call'', ''put'' and ''otm'' are allowed.']);
           end           
       end
       ttmUniq = unique(S_contracts.ttm);
       cartProdS = S_contracts.cartProd;
   end
   
   if doPriceVIX
       if ~isfield(VIX_contracts,'priceFormat')
           VIX_contracts.priceFormat = 'implied_volatility';
       end       
       q_is_used_VIX = isFieldNonEmpty(VIX_contracts,'q');
       k_is_used_VIX = isFieldNonEmpty(VIX_contracts,'k');
       K_is_used_VIX = isFieldNonEmpty(VIX_contracts,'K');
       if q_is_used_VIX + k_is_used_VIX + K_is_used_VIX ~= 1
            error(['PricingModelClass: Exactly one of the fields ''K'',''k'' and ''q'' must be ',...
                   'specified for the ''VIX_contracts'' input.']);
       end
       if q_is_used_VIX
            moneyness_VIX = VIX_contracts.q;
            if any(VIX_contracts.q <= 0 | VIX_contracts.q >= 1)
              error('PricingModelClass:GetPrices: Values of ''VIX_contracts.q'' must lie in (0,1)');
            end            
       elseif K_is_used_VIX
            moneyness_VIX = VIX_contracts.K;
       elseif k_is_used_VIX
            moneyness_VIX = VIX_contracts.k;
       end
       if isFieldEmpty(VIX_contracts,'cartProd')
           if size(moneyness_VIX,1) ~= size(VIX_contracts.ttm,1)
               VIX_contracts.cartProd = true;
           else
               VIX_contracts.cartProd = false;
           end
       end
       if ~isfield(VIX_contracts,'optionType') || isempty(VIX_contracts.optionType)
           VIX_contracts.optionType = 'otm';
       elseif islogical(VIX_contracts.optionType)
           if VIX_contracts.cartProd 
               error(['PricingModelClass:GetPrices: If the optionType ',...
                      'parameter is specified as a logical vector we must have ',...
                      'cartProd = false.']);
           elseif size(moneyness_VIX,1) ~= size(VIX_contracts.optionType,1)  || ...
                  size(moneyness_VIX,1) ~= size(VIX_contracts.ttm,1)
               error(['PricingModelClass:GetPrices: If the optionType parameter is specified ',...
                      'as a logical vector, the optionType, moneyness and maturity vectors ',...
                      'must have the same sizes.']);
           end
       else
           if ~ismember(VIX_contracts.optionType,{'call','put','otm'})
               error(['PricingModelClass:GetPrices: If the parameter ',...
                      'optionType is specified as a string only the values ',...
                      '''call'', ''put'' and ''otm'' are allowed.']);
           end           
       end
       if ~VIX_contracts.cartProd && ~(size(moneyness_VIX,1) == size(VIX_contracts.ttm,1))
           error(['PricingModelClass:GetPrices: When not using a cartesian',...
                  ' product the moneyness and time-to-expiry vectors must be', ...
                  ' of the same length.']);           
       end
       if VIX_contracts.cartProd && (any(size(unique(moneyness_VIX)) ~= size(moneyness_VIX)) ...
               || any(size(unique(VIX_contracts.ttm)) ~= size(VIX_contracts.ttm)))
            error(['PricingModelClass:GetPrices: The moneyness and expiry ',...
                   'vectors must be unique when a cartProd product is ',...
                   'requested.']);
       end
       ttmUniq = unique([ttmUniq;VIX_contracts.ttm]);
       cartProdV = VIX_contracts.cartProd;
   end
   
   if doPriceS
       if q_is_used_S
           Ns = size(S_contracts.q,1);
           S_contracts.k = [];
           S_contracts.K = [];
       elseif k_is_used_S
           Ns = size(S_contracts.k,1);
           S_contracts.q = [];
           S_contracts.K = [];
       elseif K_is_used_S
           Ns = size(S_contracts.K,1);
           S_contracts.q = [];
           S_contracts.k = [];
       end
       Ms = size(S_contracts.ttm,1);
   end
      
   if doPriceVIX
       if q_is_used_VIX
           Nv = size(VIX_contracts.q,1);
           VIX_contracts.K = [];
           VIX_contracts.k = [];
       elseif k_is_used_VIX
           Nv = size(VIX_contracts.k,1);
           VIX_contracts.q = [];
           VIX_contracts.K = [];
       elseif K_is_used_VIX
           Nv = size(VIX_contracts.K,1);
           VIX_contracts.q = [];
           VIX_contracts.k = [];
       end
       Mv = size(VIX_contracts.ttm,1);
   end
   
   % Extract pricer-settings:
   settings = obj.pricerSettings;
   
   % If we are to aggregate prices over several runs then call function recursively (can 
   % be useful to limit memory usage):
    if isa(settings,'MonteCarloPricerSettingsClass') && settings.numRepeat > 1 ...
            && ~specialRun  
        
        if useExistingNumbers
            error(['PricingModelClass:GetPrices: Having useExistingNumbers = true is not ',...
                   'compatible with obj.pricerSettings.numRepeat > 1.']);
        end
        
        % Check for any old random numbers:
        if ~isempty(settings.random)
          warning(['PricingModelClass: GetPrices: When using the ''numRepeat'' parameter of ',...
                 'the MonteCarloPricerSettingsClass there must not be any random numbers ',...
                 'stored in the object before running. Clearing existing random numbers ',...
                 'stored in obj.pricerSettings.random in order to simulate new ones.']);
          settings.random = [];
        end
        
        % Run recursively:
        nRep = settings.numRepeat;
        
        [pS,pV] = deal([]);
        S_contracts_tmp = S_contracts;
        VIX_contracts_tmp = VIX_contracts;
        for i=1:nRep
            % Get prices:
            [pS_new,pV_new,ivIdx] = obj.GetPrices(S_contracts_tmp,VIX_contracts_tmp,...
                                            'useExistingNumbers',false,'specialRun',true);
            
            if i > 1 && doPriceS && q_is_used_S && cartProdS
                % Reshape vector output to matrices:
                pS_new.p = reshape(pS_new.p,size(pS.p,1),size(pS.p,2));
                pS_new.se = reshape(pS_new.se,size(pS.se,1),size(pS.se,2));
            end
            if i > 1 && doPriceVIX && (q_is_used_VIX || k_is_used_VIX) && cartProdV
                % Reshape vector output to matrices:
                pV_new.p = reshape(pV_new.p,size(pV.p,1),size(pV.p,2));
                pV_new.se = reshape(pV_new.se,size(pV.se,1),size(pV.se,2));
                pV_new.F = reshape(pV_new.F,size(pV.F,1),size(pV.F,2));
                pV_new.se_F = reshape(pV_new.se_F,size(pV.se_F,1),size(pV.se_F,2));                
            end
            
            % Accumulate:
            if i == 1
                pS = pS_new;
                pV = pV_new;
            else
                if ~isempty(pS) && ~isempty(fieldnames(pS))
                    pS.p = pS.p + pS_new.p;
                    pS.se = pS.se + pS_new.se;
                end
                
                if ~isempty(pV) && ~isempty(fieldnames(pV))
                    if ivIdx
                        error(['PricingModelClass:GetPrices: numRepeat > 1 is currently not ',...
                               'supported for the pricing of VIX options under models where ',...
                               'implied volatility is returned from the GetPricesSub method.']);
                    end
                    % Comment: For VIX options the ATM strike is the VIX futures price which is
                    % estimated using Monte Carlo and therefore can change from iteration to
                    % iteration. We therefore need to appropriately convert prices to ensure 
                    % correctness:
                    convCall = pV.idxCall == true & pV_new.idxCall == false;
                    convPut = pV.idxCall == false & pV_new.idxCall == true;
                    if any(convCall)
                        pV_new.p(convCall) = pV_new.p(convCall) ...
                                                + (pV_new.F(convCall) - pV_new.K(convCall))...
                                                .*exp(-obj.y.Eval(pV_new.ttm(convCall)));    
                    end
                    if any(convPut)
                        pV_new.p(convPut) = pV_new.p(convPut) ...
                                                - (pV_new.F(convPut)-pV_new.K(convPut))...
                                                .*exp(-obj.y.Eval(pV_new.ttm(convPut)));                
                    end
                    
                    % Aggregate results:
                    pV.p = pV.p + pV_new.p;
                    pV.se = pV.se + pV_new.se;
                    pV.F = pV.F + pV_new.F;
                    pV.se_F = pV.se_F + pV_new.se_F;
                end
            end
            
            % If 'q' (i.e. quantiles) is used for moneyness (for S and/or VIX options), we 
            % compute the corresponding strikes from the first run and fix them from thereon. 
            % We do the same for VIX options if log-moneyness (i.e. 'k') is specified. The 
            % reason being VIX futures prices are also estimated using Monte Carlo.
            if i == 1
                if isfield(S_contracts,'q') && ~isempty(S_contracts.q)
                    if ~S_contracts.cartProd
                        S_contracts_tmp.k = pS.k;
                        S_contracts_tmp.q = [];
                    else
                        % When Cartesian product is used in this context we need to reshape
                        % to vectors to conform with input requirements:
                        S_contracts_tmp.k = pS.k(:);
                        S_contracts_tmp.ttm = pS.ttm(:);
                        S_contracts_tmp.optionType = logical(pS.idxCall(:));
                        S_contracts_tmp.q = [];  
                        S_contracts_tmp.cartProd = false;
                    end
                end
                if isfield(VIX_contracts,'q') && ~isempty(VIX_contracts.q) ...
                        || isfield(VIX_contracts,'k') && ~isempty(VIX_contracts.k)
                    if ~VIX_contracts.cartProd
                        VIX_contracts_tmp.K = pV.K;
                        VIX_contracts_tmp.k = [];
                        VIX_contracts_tmp.q = [];
                    else
                        % When Cartesian product is used in this context we need to reshape
                        % to vectors to conform with input requirements:                        
                        VIX_contracts_tmp.K = pV.K(:);
                        VIX_contracts_tmp.ttm = pV.ttm(:);
                        VIX_contracts_tmp.optionType = logical(pV.idxCall(:));
                        VIX_contracts_tmp.q = [];   
                        VIX_contracts_tmp.k = [];
                        VIX_contracts_tmp.cartProd = false;
                    end                    
                end 
            end
            disp(['PricingModelClass:GetPrices: Iteration ', num2str(i),...
                  ' out of ',num2str(nRep),' completed']);
        end
        
        if ~isempty(pS) && ~isempty(fieldnames(pS))
            pS.p = pS.p ./ nRep;
            pS.se = (pS.se ./ nRep)./sqrt(nRep);
        end
        
        if ~isempty(pV) && ~isempty(fieldnames(pV))
            pV.p = pV.p ./ nRep;
            pV.F = pV.F ./ nRep;
            pV.se = (pV.se ./ nRep)./sqrt(nRep);
            pV.se_F = (pV.se_F ./ nRep)./sqrt(nRep);
        end
        
    else
        %% Prepare random numbers:
        if isa(settings,'MonteCarloPricerSettingsClass')
           if ~useExistingNumbers
               if ~isempty(settings.random)
                  warning(['PricingModelClass:GetPrices: Clearing existing ',...
                    'random numbers in object in order to simulate new ones.']);
               end
               settings.ClearNumbers();
               obj.GenNumbersForMC(ttmUniq,doPriceS,doPriceVIX, exist('q_is_used_S','var')...
                                                                && q_is_used_S, ...
                                                                obj.pricerSettings.seed);
               clearNumbersAfterRun = true;
           elseif useExistingNumbers && ~isempty(settings.random)
               clearNumbersAfterRun = false;
           elseif useExistingNumbers && isempty(settings.random) 
             error(['PricingModelClass:GetPrices: ''useExistingNumbers''', ...
                   ' was set to true but there is no random numbers stored',...
                   ' in the object. Suggestions: Either (1) set ', ...
                   '''useExistingNumbers'' to false or (2) call the ',...
                   '''GenNumbersForMC'' method with the relevant expiries. ',...
                   'Then call the ''GetPrices'' method again.']);               
           end
        end

        %% Group the expiries appropriately:
        if isa(settings,'MonteCarloPricerSettingsClass')
            % We group the expiries according to how many steps we simulate:
            [~,~,idxn] = unique(settings.GetNumSteps(ttmUniq));
            uniq_grp = unique(idxn);

        elseif strcmpi(class(obj.pricerSettings),'FourierPricerSettingsClass')
            % Each expiry is computed separately:
             [~,~,idxn]= unique(ttmUniq); 
             uniq_grp = unique(idxn);
        end
        
        % Spread indices to each of the expiration vectors:
        if doPriceS
            [~,idxtemp] = ismembertol(S_contracts.ttm,ttmUniq);
            idxn_S = idxn(idxtemp);
        end
        if doPriceVIX
            [~,idxtemp] = ismembertol(VIX_contracts.ttm,ttmUniq);
            idxn_VIX = idxn(idxtemp);
        end        

        %% Initialize result matrices:
        if doPriceS
            % Find matrix size:
            if ~S_contracts.cartProd
                output_size = [Ns,1];
            else
                output_size = [Ns,Ms];
            end
            % Initialize:
            [pS.p,pS.idxCall,pS.F,pS.K,pS.k,pS.ttm,pS.se] = deal(NaN(output_size));
            if q_is_used_S
                if S_contracts.cartProd
                    pS.q = repmat(S_contracts.q,1,Ms);
                else
                    pS.q = S_contracts.q;
                end
            else
                pS.q = [];
            end
            pS.priceFormat = S_contracts.priceFormat;
        end    
        if doPriceVIX
            % Find matrix size:
            if ~VIX_contracts.cartProd
                output_size = [Nv,1];
            else
                output_size = [Nv,Mv];
            end
            % Initialize:
            [pV.p,pV.idxCall,pV.F,pV.K,pV.k,pV.ttm,pV.se,pV.se_F] = deal(NaN(output_size));
            if q_is_used_VIX
                if VIX_contracts.cartProd
                    pV.q = repmat(VIX_contracts.q,1,Mv);
                else
                    pV.q = VIX_contracts.q;
                end  
            else
                pV.q = [];
            end
            pV.priceFormat = VIX_contracts.priceFormat;
        end

        %% Compute prices 
        % ... and do so separately for each expiration group:
        for i=1:size(uniq_grp,1)
            % Determine whether or not to run the 'GetPricesSub' method using a cartProd 
            % product or not:
            if isa(obj.pricerSettings,'FourierPricerSettingsClass') 
                if doPriceS;cartProd_S_sub = true;end
                if doPriceVIX;cartProd_V_sub = true;end
            else
                if doPriceS;cartProd_S_sub = S_contracts.cartProd;end
                if doPriceVIX;cartProd_V_sub = VIX_contracts.cartProd;end
            end

            %% Subset maturity and moneyness inputs for this iteration:
            if doPriceS
                idxTs = idxn_S == i;
                if ~S_contracts.cartProd
                    if q_is_used_S
                        KS_sub = [];
                        kS_sub = [];
                        qS_sub = S_contracts.q(idxTs);
                    elseif k_is_used_S
                        KS_sub = [];
                        kS_sub = S_contracts.k(idxTs);
                        qS_sub = [];
                    elseif K_is_used_S
                        KS_sub = S_contracts.K(idxTs);
                        kS_sub = [];
                        qS_sub = [];
                    end
                else
                    KS_sub = S_contracts.K;
                    kS_sub = S_contracts.k;
                    qS_sub = S_contracts.q;
                end
                ttmS_sub = S_contracts.ttm(idxTs);
            else
                [KS_sub,kS_sub,qS_sub,ttmS_sub,cartProd_S_sub] = deal([]);
            end

            if doPriceVIX
                idxTv = idxn_VIX == i;
                if ~VIX_contracts.cartProd
                    if q_is_used_VIX
                        kV_sub = [];
                        KV_sub = [];
                        qV_sub = VIX_contracts.q(idxTv);
                    elseif K_is_used_VIX
                        KV_sub = VIX_contracts.K(idxTv);
                        kV_sub = [];
                        qV_sub = [];
                    elseif k_is_used_VIX
                        kV_sub = VIX_contracts.k(idxTv);
                        KV_sub = [];
                        qV_sub = [];                        
                    end
                else
                    kV_sub = VIX_contracts.k;
                    KV_sub = VIX_contracts.K;
                    qV_sub = VIX_contracts.q;
                end
                ttmV_sub = VIX_contracts.ttm(idxn_VIX == i);
            else
                [kV_sub,KV_sub,qV_sub,ttmV_sub,cartProd_V_sub] = deal([]);
            end
            
                % Run pricing algorithm for subset of contracts:
                [pS_tmp,idxCallS_tmp,seS_tmp,KS_sub_out,kS_sub_out,...
                 pV_tmp,idxCallV_tmp,seV_tmp,KV_sub_out,kV_sub_out,...
                 futV_tmp,seFutV_tmp,ivIdx] = obj.GetPricesSub(KS_sub,kS_sub,qS_sub,ttmS_sub,...
                                                               cartProd_S_sub,...
                                                               KV_sub,kV_sub,qV_sub,...
                                                               ttmV_sub,cartProd_V_sub);

            % Store results:
            if doPriceS && ~isempty(pS_tmp)
                if S_contracts.cartProd
                    pS.k(:,idxTs') = kS_sub_out;
                    pS.K(:,idxTs') = KS_sub_out;
                    pS.ttm(:,idxTs') = repmat(ttmS_sub',size(kS_sub_out,1),1);
                    pS.p(:,idxTs') = pS_tmp;
                    pS.idxCall(:,idxTs') = idxCallS_tmp;
                    pS.se(:,idxTs')= seS_tmp;
                else
                    pS.k(idxTs) = kS_sub_out;
                    pS.K(idxTs) = KS_sub_out;
                    pS.ttm(idxTs) = ttmS_sub;
                    pS.p(idxTs) = pS_tmp;
                    pS.idxCall(idxTs) = idxCallS_tmp;
                    pS.se(idxTs) = seS_tmp;
                end
            end

            if doPriceVIX && ~isempty(pV_tmp)
                if VIX_contracts.cartProd
                    pV.K(:,idxTv') = KV_sub_out; 
                    pV.k(:,idxTv') = kV_sub_out; 
                    pV.ttm(:,idxTv') = repmat(ttmV_sub',size(KV_sub_out,1),1);
                    pV.p(:,idxTv') = pV_tmp;
                    pV.idxCall(:,idxTv') = idxCallV_tmp;
                    pV.se(:,idxTv') = seV_tmp;
                    pV.F(:,idxTv') = futV_tmp;
                    pV.se_F(:,idxTv') = seFutV_tmp;
                else
                    pV.K(idxTv) = KV_sub_out;
                    pV.k(idxTv) = kV_sub_out;
                    pV.ttm(idxTv) = ttmV_sub;
                    pV.p(idxTv) = pV_tmp;
                    pV.idxCall(idxTv) = idxCallV_tmp;
                    pV.se(idxTv) = seV_tmp;
                    pV.F(idxTv) = futV_tmp;
                    pV.se_F(idxTv) = seFutV_tmp;
                end
            end

        end

        % Clear random numbers (if needed):
        if isa(obj.pricerSettings,'MonteCarloPricerSettingsClass') && clearNumbersAfterRun
            settings.ClearNumbers();
        end
        
        if specialRun
            return;
        end
    end

    % Transform price output:
    if doPriceS
        pS.priceFormat = S_contracts.priceFormat;
        isNegative = pS.p <= 0;
        pS.p(isNegative) = 0;
        pS.se(isNegative) = NaN;
        pS = obj.TransformPriceOutput(pS,S_contracts,false,ivIdx);
    end
    
    if doPriceVIX
        pV.priceFormat = VIX_contracts.priceFormat;
        isNegative = pV.p <= 0;
        pV.p(isNegative) = 0;
        pV.se(isNegative) = NaN;        
        pV = obj.TransformPriceOutput(pV,VIX_contracts,true,ivIdx);
    end

end
function pOpt = TransformPriceOutput(obj,pOpt,contracts,vix,ivIdx)
%
%       Transforms prices from one format to another. Only intended for use by the 'GetPrices' 
%       method; not by an end-user.
%
% ------------------------------------------------------------------------------------------------
% Parameters
% ------------------------------------------------------------------------------------------------
%   pOpt:      [1x1 struct]  The 'pS' or 'pVIX' struct from the 'GetPrices' method.
%   contracts: [1x1 struct]  The 'S_contracts' or 'VIX_contracts' input from the 'GetPrices' 
%                            method.
%   vix:       [1x1 logical] True if input is for VIX options.
%   ivIdx:     [1x1 logical] True if price format of pOpt is 'implied_volatility'. Set to false
%                            if price format is 'monetary'.
%
% ------------------------------------------------------------------------------------------------
% Output
% ------------------------------------------------------------------------------------------------
%   pOpt:      [1x1 struct] The 'pOpt' input variable transformed.
%
% ------------------------------------------------------------------------------------------------

    % Compute additional moneyness formats:
    y = obj.y.Eval(pOpt.ttm);
    if ~vix
        q = obj.q.Eval(pOpt.ttm);
        pOpt.F = obj.s0.*exp((y-q).*pOpt.ttm);
        pOpt.K = pOpt.F.*exp(pOpt.k);
    else
        pOpt.k = log(pOpt.K./pOpt.F);
    end
    
    idxCall = pOpt.idxCall;
    
    % Compute and/or transform price format and/or price type:
    if strcmpi(pOpt.priceFormat,'implied_volatility') && ~ivIdx
        % Compute implied volatilities
        if ~vix
            pOpt.p = blsimpv_with_negative_rates(obj.s0,pOpt.K,y,pOpt.ttm,q,pOpt.p,pOpt.idxCall);
        else
            pOpt.p = blsimpv_with_negative_rates(pOpt.F,pOpt.K,0,pOpt.ttm,...
                                                 0,pOpt.p.*exp(y.*pOpt.ttm),pOpt.idxCall);
        end
        
        % Convert standard errors from prices to implied volatilities with delta method:
        if ~vix
            bs_vega = BSGreek('vega',[],obj.s0,pOpt.K,y,pOpt.ttm,pOpt.p,q);
        else
            bs_vega = BSGreek('vega',[],pOpt.F,pOpt.K,0,pOpt.ttm,pOpt.p.*exp(y.*pOpt.ttm),0);
        end
        pOpt.se = pOpt.se./bs_vega;
        
    elseif strcmpi(pOpt.priceFormat,'monetary') && ivIdx
        % Convert from implied volatility to prices:
        if ~vix
            pOpt.p(idxCall) = bscall(obj.s0,pOpt.K(idxCall),y(idxCall),pOpt.ttm(idxCall),...
                                     pOpt.p(idxCall).^2,q(idxCall));
            pOpt.p(~idxCall) = bsput(obj.s0,pOpt.K(~idxCall),y(~idxCall),pOpt.ttm(~idxCall),...
                                     pOpt.p(~idxCall).^2,q(~idxCall));
        else
            error(['PricingModelClass:TransformPriceOutput: Conversion from implied ',...
                   'volatility is currently not supported for VIX options.']);
        end
        
    elseif strcmpi(pOpt.priceFormat,'monetary') && ~ivIdx
        % Convert prices from puts to calls and vice versa:
        optionType = contracts.optionType;
        if ~islogical(optionType)
            switch optionType
                case 'otm'
                    convToCall = pOpt.k >= 0;
                case 'call'
                    convToCall = true(size(pOpt.k));
                case 'put'
                    convToCall = false(size(pOpt.k));
            end
        else
            convToCall = optionType;
        end

        % Is call but should be put:
        convPut = ~convToCall & idxCall;

        % Is put but should be call:
        convCall = convToCall & ~idxCall;

        % Convert:
        if ~vix
            pOpt.p(convPut) = pOpt.p(convPut) - obj.s0.*exp(-q(convPut).*pOpt.ttm(convPut)) ...
                              + pOpt.K(convPut).*exp(-y(convPut).*pOpt.ttm(convPut));
            pOpt.p(convCall) = pOpt.p(convCall) + obj.s0.*exp(-q(convCall).*pOpt.ttm(convCall)) ...
                               - pOpt.K(convCall).*exp(-y(convCall).*pOpt.ttm(convCall));
        else
            pOpt.p(convPut) = pOpt.p(convPut) - (pOpt.F(convPut) - pOpt.K(convPut))...
                                                   .*exp(-obj.y.Eval(pOpt.ttm(convPut)));
            pOpt.p(convCall) = pOpt.p(convCall) + (pOpt.F(convCall) - pOpt.K(convCall))...
                                                   .*exp(-obj.y.Eval(pOpt.ttm(convCall)));
        end
    end
    
    % Reset some unnecessary fields:
    if strcmpi(pOpt.priceFormat,'implied_volatility')
        pOpt.optionType = [];
        pOpt.idxCall = [];
    end

end
[p, idxCall, se] = MonteCarloEstimation(~,k,t,cartProd,paths,optType,condMC,cv,...
                                        antithetic,prec,F,y,rho)        
function aStar = GetOptimalFourierAlpha(obj,K,T)
%
%    Returns an optimal dampening parameter as suggested in (Lord and Kahl, 2006). We assume 
%    zero interest rates and dividends, and an initial asset price of 1.
%
%    Remarks:
%       - We use the rule of thumb suggested by (Lord and Kahl, 2006) where we for 
%         out-of-the-money calls only look for dampening parameters (alpha's) above 0 
%         and for out-of-the-money puts only look for values below -1.
%
% -------------------------------------------------------------------------------------------------
%  Parameters                                           
% -------------------------------------------------------------------------------------------------
%   K:               [Nx1 real]       Strikes.
%   T:               [1x1 real]       Expiry.
%
% -------------------------------------------------------------------------------------------------
%  Output                                               
% -------------------------------------------------------------------------------------------------
%   aStar:           [Nx1 real]       Optimal dampening parameter.
%
% -------------------------------------------------------------------------------------------------
%  References
% -------------------------------------------------------------------------------------------------
%   - Lord, R. and Kahl, C., Optimal Fourier Inversion in Semi-Analytical Option Pricing. 
%     2006, Tinbergen Institute Discussion Paper No. 2006-066/2.
%
% -------------------------------------------------------------------------------------------------

    %% Find valid range of alpha's to optimize over:
    aStar = NaN(size(K));
    idxCall = K>=1;
    use_default = false;
    [aMin,aMax] = deal(NaN);
    try 
        [aMin, aMax] = obj.GetMomentExplosion(T);
    catch
        use_default = true;
    end
    if isnan(aMin) || isnan(aMax) || abs(aMin)==Inf || abs(aMax)==Inf
        use_default = true;
    end
    if use_default
       warning(['PricingModelClass:GetOptimalFourierAlpha: Something went',...
                ' wrong while computing the moment bounds. Defaulting to a',...
                ' dampening parameter of 1/2 for out-of-the-money calls and',...
                ' -1/2 for out-of-the-money puts.']);
        % Probably assumptions were not met, we set a default:
        aStar(idxCall) = 0.5;
        aStar(~idxCall) = -0.5;
        return;
    end
    alphaMin = aMin - 1;
    alphaMax = aMax - 1;

    % Define functions:
    if isa(obj,'rHestonClass')
        xiVals = obj.CharacteristicFun(T,1,'only_final_timepoint',[],true);
        phi = @(u)(obj.CharacteristicFun(T,u,'only_final_timepoint',xiVals));
    else
        phi = @(u)(obj.CharacteristicFun(T,u));
    end

    if isa(obj,'rHestonClass') ...
            && strcmpi(obj.charFunSettings.method,'RationalApprox')
    % In this case we need to limit the alpha bounds further, see
    % (Gatheral and Radoicic, 2019).
        alphaMin = max([alphaMin;-1]);
        alphaMax = min([(1 - obj.rho.^2).^(-1) - 1;alphaMax]);
    end

    %% Define objective function:
    psi = @(v,alpha) ( real(phi( v - (alpha + 1)*1i  ) ...
                        ./ ( alpha.^2 + alpha - v.^2 + 1i*(2*alpha + 1).*v) ) );
    PSI = @(alpha,K)( -alpha*log(K).' + 0.5*log( psi(0,alpha).^2 ) );  

    %% Find optimal alphas:
    
    % Optimization settings:
    eps = 0.0001;
    options = optimset('Display','off','MaxIter',10^6,'MaxFunEvals',10^6);
    
    % Loop over the strikes:
    for i=1:size(K,1)
        
        % Adjust the optimization range according to the suggested
        % 'rule-of-thumb' and set the initial guess:
        if K(i) >= 1
            lb = 0 + eps;
            ub = alphaMax - eps;
            alpha0 = mean([0,alphaMax]);
        else
            lb = alphaMin + eps;
            ub = -1 - eps;
            alpha0 = mean([-1,alphaMin]);
        end
        
        % Readjust the initial guess again if the optimization function
        % cannot properly evaluate the initial guess:
        if isnan(PSI(alpha0,K(i))) || abs(PSI(alpha0,K(i))) == Inf
            % Try a different choice:
            if K(i) >= 1
                alpha0 = alphaMax ./ 10;
            else
                alpha0 = -1 + (alphaMin - (-1)) ./ 10;
            end

            % If there is still a problem we use the default value as an 
            % initial guess:
            if isnan(PSI(alpha0,K(i))) || abs(PSI(alpha0,K(i))) == Inf
                alpha0 = (K(i) >= 1)*0.5 - (K(i) < 1)*0.5;
                if alpha0 < 0
                    lb = -1 + eps;
                    ub = 0 - eps;
                end
            end
        end
        
        % Solve the optimization problem:
        large_number = 10^(20);
        [aOpt,~,flag] = fminsearch(@(x)(PSI(x,K(i))*(1 - (x < lb | x > ub))...
                              + (x < lb | x > ub)*large_number),alpha0,options);

        if flag ~= 1
            aStar(i) = (K(i) >= 1)*0.5 - (K(i) < 1)*0.5;
            warning(['PricingModelClass:GetOptimalFourierAlpha: Something ',...
                     'went wrong when searching for an optimal alpha ',...
                     'for the strike ', num2str(K(i)), '. Defaulting to a ',...
                     'dampening parameter of ', num2str(aStar(i))]);
        else
            aStar(i) = aOpt;
        end
        
    end

end
function [lBound, uBound] = GetMomentExplosion(obj,T,maxIter)
% 
%    Returns the lower and upper moment explosions at a given future time point T. That is, 
%    the function will return a value 'lBound' s.t. E[S(T)^(mu)] is infinite for all 
%    mu <= lBound and also return a value 'uBound' s.t. E[S(T)^(mu)] is infinite for all 
%    mu >= uBound.
%
%    The pricing model must implement the following methods:
%
%       o MomentExplosionInitialGuess:  A function with no arguments that should return valid 
%                                       initial guesses for the moment explosions.
%
%       o GetExplosionTime:             A function taking a moment as input and then returning 
%                                       the moment explosion time. This method should be 
%                                       vectorized.
%
% -------------------------------------------------------------------------------------------------
%  Parameters                                           
% -------------------------------------------------------------------------------------------------
%   T:               [Nx1 real]       Future time points to consider the distribution of S(T) at.
%   maxIter:         [1x1 integer]    Maximum number of iterations. Is optional. Default is 10^5.
%
% -------------------------------------------------------------------------------------------------
%  Outputs
% -------------------------------------------------------------------------------------------------
%   lBound:          [Nx1 real]       Lower bounds for moment explosions.
%   uBound:          [Nx1 real]       Upper bounds for moment explosions.
%
% -------------------------------------------------------------------------------------------------

   if ~any(strcmpi(methods(obj),'MomentExplosionInitialGuess'))
       error(['PricingModelClass:GetMomentExplosion: Could not find method',...
              ' ''MomentExplosionInitialGuess''.']);
   end
   if ~any(strcmpi(methods(obj),'GetExplosionTime'))
       error(['PricingModelClass:GetMomentExplosion: Could not find method',...
              ' ''GetExplosionTime''.']);
   end   

   % Set the error tolerance:
   eps = 0.0001;
   
   % Maximum number of steps:
   if ~exist('maxSteps','var') || isempty(maxIter)
       maxIter = 10^5;
   end
   
   [lBound,uBound] = deal(NaN(size(T)));
   for i=1:size(T,1)
       % Get initial guesses
       [omegaLow, omegaHigh] = obj.MomentExplosionInitialGuess();
        omegaLow = omegaLow + eps;
        omegaHigh = omegaHigh - eps;
       
       % Set the initial step size:
       dT = 1;
       
       % Reset the number of steps taken:
       nIter = 0;

       % Lower bound:
       errFun = @(omega)(T(i)-obj.GetExplosionTime(omega));
       sgn = sign(errFun(omegaLow));
       errLatest = errFun(omegaLow);
       while abs(errLatest) > eps
           % Take new step:
           omegaLow = omegaLow + sgn*dT;
           nIter = nIter + 1;
           if nIter > maxIter
               error(['PricingModelClass:GetMomentExplosion: ',...
                      'Maximum number of iterations reached.']);
           end
           while sign(errFun(omegaLow)) ~= sign(errLatest)
               % We went too far, so step back:
               omegaLow = omegaLow - sgn*dT;
               % Adjust dT
               dT = dT/10;
               % Try smaller step:
               omegaLow = omegaLow + sgn*dT;
               % Check if within error tolerance:
               if abs(errFun(omegaLow)) < eps
                   break;
               end
           end
           errLatest = errFun(omegaLow);
       end
       lBound(i) = omegaLow;

       % Reset step size:
       dT = 1;
       
       % Reset the number of steps taken:
       nIter = 0;       

       % Upper bound:
       errFun = @(omega)(T(i)-obj.GetExplosionTime(omega));
       sgn = -sign(errFun(omegaHigh));
       errLatest = errFun(omegaHigh);
       while abs(errLatest) > eps
           % Take new step:
           omegaHigh = omegaHigh + sgn*dT;
           nIter = nIter + 1;
           if nIter > maxIter
               error(['PricingModelClass:GetMomentExplosion: ',...
                      'Maximum number of iterations reached.']);
           end           
           while sign(errFun(omegaHigh)) ~= sign(errLatest)
               % We went too far, so step back:
               omegaHigh = omegaHigh - sgn*dT;
               % Adjust dT
               dT = dT/10;
               % Try smaller step:
               omegaHigh = omegaHigh + sgn*dT;
               % Check if within error tolerance:
               if abs(errFun(omegaHigh)) < eps
                   break;
               end
           end
           errLatest = errFun(omegaHigh);
       end    
       uBound(i) = omegaHigh;
   end

end
function result = Calibrate(obj,S_contracts,VIX_contracts,parNames,par0,lb,ub,varargin)    
%
%    Calibrates to European calls and puts by minimising the following error measure
%
%         wRMSE :=  sqrt( sum_{i=1}^n w(i)*err(i)^2 )
%
%    where 'w' are weights and 'err' either is the error between model implied volatilities and
%    mid implied volatilities or the deviation of model values from the bid-ask spread.
%
%    Weights are automatically normalised so sum(w) = 1.
%
%    Object properties will be changed during the optimisation; if succesful they are set to 
%    their optimal values. 
% 
% -----------------------------------------------------------------------------------------------
%  Main parameters  
% -----------------------------------------------------------------------------------------------
%   S_contracts:     [1 x 1 struct]   Information on S(t) options. The following fields apply:
%       o k:         [Ns x 1 real]    Log-moneyness values.
%       o ttm:       [Ns x 1 real]    Maturities.
%       o mid:       [Ns x 1 real]    Implied volatility (mid). 
%       o bid:       [Ns x 1 real]    Implied volatility (bid).
%       o ask:       [Ns x 1 real]    Implied volatility (ask).
%       o w:         [Ns x 1 real]    Weights. Must be positive. Defaults to a uniform vector 
%                                     if empty or unspecified. 
%
%   VIX_contracts:   [1 x 1 struct]   Information on VIX options. The following fields apply:
%       o K:         |Nv x 1 real]    Strikes.
%       o ttm:       [Nv x 1 real]    Maturities.
%       o mid:       [Nv x 1 real]    Implied volatility (mid).
%       o bid:       [Nv x 1 real]    Implied volatility (bid).
%       o ask:       [Nv x 1 real]    Implied volatility (ask).
%       o w:         [Nv x 1 real]    Weights. Must be positive. Defaults to a uniform vector 
%                                     if empty or unspecified. 
% 
%   parNames:        [1xN cell]       Object property (parameter) names. Must either refer to   
%                                     scalars or CurveClass objects.
%
%   par0:            [1xN cell]       Initial guess. Each element should be a column vector with  
%                                     length compatible with the corresponding object property.
%
%   lb:              [1xN cell]       Lower bounds. Specify similarly to 'par0'.
%
%   ub:              [1xN cell]       Upper bounds. Specify similarly to 'par0'.
%
% -----------------------------------------------------------------------------------------------
%  Optional parameters  
% -----------------------------------------------------------------------------------------------
% * Specify as name-value pairs.
%
%   options:         [1x1 Lsqnonlin]  Options used by the 'Lsqnonlin' optimisation algorithm.
%                                     See code for default.
%
%   error_measure:   [1x1 string]     Options are 'dist_to_mid' and 'dist_to_spread'.
%
%   dispIter:        [1x1 logical]    Shows parameter values in each iteration. Default is false.
%
%   plotFit:         [1x1 logical]    Plots fit for each iteration. Default is false.
%
%   nS:              [1x1 integer]    Number of plots to show for options on S when plotFit = true.
%                                     Default is 5.
%
%   nV:              [1x1 integer]    Number of plots to show for options on V when plotFit = true.
%                                     Default is 5.
%
%   parEqual:        [1xM cell]       Parameters which should be held equal to other parameters. 
%                                     Each element is a [1x2 cell] of strings, say, {'param1',
%                                     'param2'}. After updating the parameter values associated
%                                     with 'parNames' we then set obj.('param1') = obj.('param2') 
%                                     for all such pairs. 
% 
% ------- Parameter group ---------
%   fSet:            [1xL cell]       Each element is a [1x1 function]. In each iteration, after 
%                                     updating parameter values associated with 'parNames' and 
%                                     'parEqual', we for i=1,...,L, set 
%                                          obj.(fOutPropNames{i}{j}) = fSet(obj,par_i1,...,par_im)
%                                     for j=1,...,length(fOutPropNames{i}), with par_i1,...,par_im 
%                                     being the parameter values associated with fInParNames{i}.
%
%   fInParNames:     [1xL cell]       Each element is a cell array of strings being a subset of
%                                     those in 'parNames'.
%
%   fOutPropNames:   [1xL cell]       Each element is a cell array consisting of strings referring
%                                     to properties of 'obj'.
%
% -----------------------------------------------------------------------------------------------
%  Outputs
% -----------------------------------------------------------------------------------------------
%   result:          [1x1 struct]     All relevant information regarding the calibration.
%
% -------------------------------------------------------------------------------------------------

    % Parse variable inputs:
    p = inputParser;
    addParameter(p,'options',optimoptions('lsqnonlin','Algorithm','trust-region-reflective',...
                                          'Display','iter','DiffMinChange',0.001));
    addParameter(p,'error_measure','dist_to_mid');
    addParameter(p,'dispIter',false);
    addParameter(p,'plotFit',false);
    addParameter(p,'nS',5);
    addParameter(p,'nV',5);
    addParameter(p,'parEqual',[]);
    addParameter(p,'fSet',[]);
    addParameter(p,'fInParNames',[]);
    addParameter(p,'fOutPropNames',[]);
    parse(p,varargin{:});
    v2struct(p.Results);
    
    % Collect weights and prices:
    [mid,bid,ask,w,moneyness,ttm,idxS] = deal([]);
    if ~isempty(S_contracts)
        ttm = [ttm;S_contracts.ttm];
        moneyness = [moneyness;S_contracts.k];
        mid = [mid;S_contracts.mid];
        if isfield(S_contracts,'bid') && ~isempty(S_contracts.bid)
            bid = [bid;S_contracts.bid];
        end
        if isfield(S_contracts,'ask') && ~isempty(S_contracts.ask)
            ask = [ask;S_contracts.ask];
        end        
        if isfield(S_contracts,'w') && ~isempty(S_contracts.w)
            w = [w;S_contracts.w];
        else
            w = [w;ones(size(S_contracts.mid))];
        end
        idxS = [idxS;ones(size(S_contracts.mid))];
        S_contracts = rmfield(S_contracts,'mid'); 
    end
    if ~isempty(VIX_contracts)
        ttm = [ttm;VIX_contracts.ttm];
        moneyness = [moneyness;VIX_contracts.K];
        mid = [mid;VIX_contracts.mid];
        if isfield(VIX_contracts,'bid') && ~isempty(VIX_contracts.bid)
            bid = [bid;VIX_contracts.bid];
        end
        if isfield(VIX_contracts,'ask') && ~isempty(VIX_contracts.ask)
            ask = [ask;VIX_contracts.ask];
        end
        if isfield(VIX_contracts,'w') && ~isempty(VIX_contracts.w)
            w = [w;VIX_contracts.w];
        else
            w = [w;ones(size(VIX_contracts.mid))];
        end
        idxS = [idxS;zeros(size(VIX_contracts.mid))];
        VIX_contracts = rmfield(VIX_contracts,'mid');
    end
    idxS = logical(idxS);
    
    % Validate and rescale weights:
    if any(w < 0)
        error('PricingModelClass:Calibrate: Weights must be non-negative.');
    end
    w = w./sum(w);

    % Validate other inputs:
    if ~isSameSize(parNames,par0,lb,ub)
        error(['PricingModelClass:Calibrate: Inputs ''parNames'', ''par0'', ''lb'' and ''ub''',...
               ' must have the same sizes.']);
    end    
    
    if ~isSameSize(fSet,fInParNames,fOutPropNames)
        error(['PricingModelClass:Calibrate: Inputs ''fSet'', ''fInParNames''',...
               ' and ''fOutPropNames'' must have the same sizes.']);
    end
    
    if ~isempty(S_contracts) 
        if isfield(S_contracts,'w')
            S_contracts = rmfield(S_contracts,'w');
        end
        if isfield(S_contracts,'bid')
            S_contracts = rmfield(S_contracts,'bid');
        end
        if isfield(S_contracts,'ask')
            S_contracts = rmfield(S_contracts,'ask');
        end        
        if ~isStructFieldsSameSize(S_contracts)
           error('PricingModelClass:Calibrate: ''S_contracts'' fields must be of the same size.');
        end
        S_contracts.priceFormat = 'implied_volatility';
    end
    
    if ~isempty(VIX_contracts) 
        if isfield(VIX_contracts,'w')
            VIX_contracts = rmfield(VIX_contracts,'w');
        end
        if isfield(VIX_contracts,'bid')
            VIX_contracts = rmfield(VIX_contracts,'bid');
        end
        if isfield(VIX_contracts,'ask')
            VIX_contracts = rmfield(VIX_contracts,'ask');
        end                
        if ~isStructFieldsSameSize(VIX_contracts)
           error('PricingModelClass:Calibrate: ''VIX_contracts'' fields must be of the same size.');
        end
        VIX_contracts.priceFormat = 'implied_volatility';
    end

    % Unpack initial guess and bounds and construct a key between parameter values 
    % and parameter names: 
    [par0_adj,lb_adj,ub_adj,parIdx] = deal([]);
    for i=1:size(parNames,2)
        par0_adj = [par0_adj;par0{i}];
        lb_adj = [lb_adj;lb{i}];
        ub_adj = [ub_adj;ub{i}];
        
        % Construct key indices:
        if i==1
            newIdx = 1;
        else
            newIdx = max(parIdx)+1;
        end
        parIdx = [parIdx;newIdx*ones(size(par0{i}))];
        
    end
    
    % Compute expanded parNames vector:
    parNamesExpanded = {};
    for i=1:size(parNames,2)
        if isprop(obj,parNames{i}) && strcmpi(class(obj.(parNames{i})),'CurveClass')
            for j=1:sum(parIdx==i)
                parNamesExpanded = {parNamesExpanded{:},[parNames{i},'(',num2str(j),')']};
            end
        else
            parNamesExpanded = {parNamesExpanded{:},parNames{i}};
        end
    end 
	
    % Compute expiries to plot (if so desired):
    [ttm_S_plot,ttm_V_plot] = deal([]);
    if plotFit
        if ~isempty(S_contracts)
            uniqT = unique(S_contracts.ttm);
            nS = min(nS,size(uniqT,1));
            dT = floor(size(uniqT,1)/nS);
            ttm_S_plot = uniqT(1:dT:end);
            ttm_S_plot = ttm_S_plot(1:nS);
            ttm_S_plot(end) = uniqT(end);
        end
        if ~isempty(VIX_contracts)
            uniqT = unique(VIX_contracts.ttm);
            nV = min(nV,size(uniqT,1));
            dT = floor(size(uniqT,1)/nV);
            ttm_V_plot = uniqT(1:dT:end);
            ttm_V_plot = ttm_V_plot(1:nV);
            ttm_V_plot(end) = uniqT(end);
        end        
    end
        
	% Fix random numbers:
	obj.GenNumbersForMC(unique(ttm),~isempty(S_contracts),~isempty(VIX_contracts),...
                        exist('S_contracts','var') && isfield(S_contracts,'q') ...
                        && ~isempty(S_contracts.q),obj.pricerSettings.seed);

    % Define pricing function: 
    p_fun = @()(obj.GetPrices(S_contracts,VIX_contracts,'useExistingNumbers',true));
    p_fun_unpacked = @()(obj.UnpackPricesForCalibration(p_fun));
    iv_model = @(par)(obj.ModifyAndGetPrices(p_fun_unpacked,parNames,par,parIdx,...
                                             parEqual,fSet,fInParNames,fOutPropNames));

    % Add output function to optimizer settings:
    obj.pricerSettings.optim = [];
    options.OutputFcn = {@(x,optimValues,state)(obj.CalibrationOutputFcn(dispIter,x,...
                                                                        parNamesExpanded,...
                                                                        optimValues,state,...
                                                                        w,idxS));...
                         @(x,optimValues,state)(obj.CalibrationPlotFcn(plotFit,state,moneyness,...
                                                                    ttm,ttm_S_plot,ttm_V_plot,...
                                                                    bid,ask,idxS,w))};
                      
    % Define error function:                                                                    
    if strcmpi(error_measure,'dist_to_mid')                                         
        err_fun = @(par)(sqrt(w).*obj.dist_to_mid(mid,iv_model(par)));
    elseif strcmpi(error_measure,'dist_to_spread')
        err_fun = @(par)(sqrt(w).*obj.dist_to_spread(bid,ask,iv_model(par)));
    else
        error('PricingModelClass:Calibrate: Invalid error measure.');
    end
    
    % Calibrate model:
    tic;
    [par_opt_vec,~,~,exitflag,output] = lsqnonlin(err_fun,par0_adj,lb_adj,ub_adj,options);
    compTime = toc;         
    
    % Get results from temporary model object (and reset):
    result = obj.pricerSettings.optim;
    obj.pricerSettings.optim = [];    
    
    % Add final pieces of information to 'result' object:
    result.exitflag = exitflag;
    result.iter = output.iterations;
    result.FCount = output.funcCount;
    result.compTime = compTime;
    
    result = rmfield(result,'iv_model');
    result = rmfield(result,'error');
    
    % Compute final error:
    ivFit = iv_model(par_opt_vec);
    if strcmpi(error_measure,'dist_to_spread')
        unw_err_sqr = obj.dist_to_spread(bid,ask,ivFit).^2;     
    elseif strcmpi(error_measure,'dist_to_mid')
        unw_err_sqr = abs(mid-ivFit).^2;
    end
    result.iv_model = ivFit;
    result.mid = mid;
    result.bid = bid;
    result.ask = ask;
    result.wRMSE = sqrt(sum(w.*unw_err_sqr));
    result.wRMSE_S = sqrt(sum((w(idxS)./sum(w(idxS))).*unw_err_sqr(idxS)));
    result.wRMSE_VIX = sqrt(sum((w(~idxS)./sum(w(~idxS))).*unw_err_sqr(~idxS)));    
    
    % Unpack parameters to a cell array:
    parOpt = cell(size(parNames));
    for i=1:size(parNames,2)
        parOpt{i} = par_opt_vec(parIdx==i);
    end
    result.par = parOpt';
    
    % Clear random numbers:
    obj.pricerSettings.random = [];    

end
function err = dist_to_mid(obj,mid,ivModel)
%
%    Computes calibration errors as distance to mid implied volatility.
%

    err = abs(mid - ivModel);
    obj.pricerSettings.optim.iv_model = ivModel;
    obj.pricerSettings.optim.error = err;
end
function err = dist_to_spread(obj,bid,ask,ivModel)
%
%    Computes calibration errors as distance to bid-ask spread (in implied volatility).
%
    err = zeros(size(bid));
    idxAbove = ivModel > ask;
    idxBelow = ivModel < bid;
    err(idxAbove) = ivModel(idxAbove) - ask(idxAbove);
    err(idxBelow) = bid(idxBelow) - ivModel(idxBelow);
    obj.pricerSettings.optim.iv_model = ivModel;
    obj.pricerSettings.optim.error = err;
end
function p = UnpackPricesForCalibration(~,p_fun)
%
%    Auxiliary function used by the 'Calibrate' method to unpack pricing data from 
%    a call to the 'GetPrices' method.
%
% -----------------------------------------------------------------------------------------------
%  Parameter
% -----------------------------------------------------------------------------------------------
%   p_fun:           [1x1 function]   Pricing function with no arguments.
%
% -----------------------------------------------------------------------------------------------
%  Output
% -----------------------------------------------------------------------------------------------
%   p:               [Nx1 real]       Unpacked vector of prices.
%
% -----------------------------------------------------------------------------------------------

    [pS,pVIX] = p_fun();
    
    p = [];
    if length(fieldnames(pS)) > 0
        p = [p;pS.p];
    end
    if length(fieldnames(pVIX)) > 0
        p = [p;pVIX.p];
    end

end
function stop = CalibrationPlotFcn(obj,plotFit,state,moneyness,ttm,ttm_S_plot,ttm_V_plot,...
                                   bid,ask,idxS,w)
%
%    'Plot function' for calibration optimizer.
%
    
    stop = false;
    if ~plotFit
        return;
    end
    
    nS = size(ttm_S_plot,1);
    nV = size(ttm_V_plot,1);
    ttmPlot = [ttm_S_plot;ttm_V_plot];
    
    if strcmpi(state,'init')
        % Create figures:
        for i=1:(nS+nV)
            fig = figure;
            fig.Name = ['calibration_',num2str(i)];
        end
        autoArrangeFigures();
        
    elseif strcmpi(state,'iter')
        % Update figures:
        figs = findobj('Type', 'figure');
        w_err_sqr = w.*obj.pricerSettings.optim.error.^2;
        for i=1:(nS+nV)
            for j=1:size(figs,1)
                if strcmpi(figs(j).Name,['calibration_',num2str(i)])
                    set(0, 'currentfigure', figs(j)); 

                    cla reset;
                    
                    if i <= nS
                        % Option on S(t):
                        idx = ismember(ttm,ttmPlot(i)) & idxS;
                        moneyness_label = 'Log-moneyness';
                    else
                        % Option on VIX(t):
                        idx = ismember(ttm,ttmPlot(i)) & ~idxS;
                        moneyness_label = 'Strike';
                    end
                    
                    hold on;plot(moneyness(idx),bid(idx),'.','color','red');
                    hold on;plot(moneyness(idx),ask(idx),'.','color','blue');
                    hold on;plot(moneyness(idx),obj.pricerSettings.optim.iv_model(idx),...
                                 '-','color','green');
                    xlabel(moneyness_label,'interpreter','latex');
                    ylabel('Implied volatility','interpreter','latex');
                    perc_of_error = sum(w_err_sqr(idx))./sum(w_err_sqr);
                    title({['T = ',num2str(ttmPlot(i))],[' (perc. of error = ',...
                            num2str(round(perc_of_error,2)),')']},'interpreter','latex');
                    
                    yyaxis right;
                    
                    hold on;bar(moneyness(idx),w_err_sqr(idx));
                    ylabel('(weighted) squared error');
                    ylim([0,5*max(w_err_sqr)]);
                    
                    if size(obj.pricerSettings.optim.wRMSE_history,1) == 1
                        allAx = findall(figs(j),'type','axes');
                        allAx.Position(2) = allAx.Position(2) - 0.025;
                    end
                end
            end
        end
    end
    drawnow;
    
end
function stop = CalibrationOutputFcn(obj,dispPar,par,parNames,optimValues,state,w,idxS)
%
%    'Output function' for calibration optimizer.
%

    stop = false; 
    res_obj = obj.pricerSettings.optim;
    if strcmpi(state,'init')
        % Initialize result object:
        res_obj.iter_history = [];
        res_obj.parNames = parNames';
        res_obj.par_history = [];        
        res_obj.wRMSE_history = [];        
        res_obj.wRMSE_S_history = [];        
        res_obj.wRMSE_VIX_history = []; 
        res_obj.FCount_history = [];
        res_obj.FOC_history = [];        
        res_obj.IV_history = [];                
        res_obj.w_errors_history = [];
    elseif strcmpi(state,'iter')
        res_obj.iter_history = [res_obj.iter_history;optimValues.iteration];
        res_obj.par_history = [res_obj.par_history,par];
        res_obj.wRMSE_history = [res_obj.wRMSE_history;sqrt(optimValues.resnorm)];  
        if any(idxS)
            res_obj.wRMSE_S_history = [res_obj.wRMSE_S_history;...
                                    sqrt(sum((w(idxS)./sum(w(idxS))).*res_obj.error(idxS).^2))];
        end
        if any(~idxS)
            res_obj.wRMSE_VIX_history = [res_obj.wRMSE_VIX_history;...
                                    sqrt(sum((w(~idxS)./sum(w(~idxS))).*res_obj.error(~idxS).^2))];
        end
        res_obj.FCount_history = [res_obj.FCount_history;optimValues.funccount];
        res_obj.FOC_history = [res_obj.FOC_history;optimValues.firstorderopt];
        res_obj.IV_history = [res_obj.IV_history,res_obj.iv_model];
        res_obj.w_errors_history = [res_obj.w_errors_history,w.*res_obj.error.^2];
        
        % Display parameters in command window:
        if dispPar
            space = '   ';
            
            % Set format-specification:
            formatspec = ['Iter = %i',space,'F-count = %i',space,...
                          'FOC = %.2e',space,'wRMSE = %i bps'];
            for i=1:size(parNames,2)
                formatspec = [formatspec,space,parNames{i},' = %8.4f'];
            end
            formatspec = [formatspec,' \n'];            
            
            % Print:
            fprintf(formatspec,[optimValues.iteration,optimValues.funccount,...
                                optimValues.firstorderopt,...
                                round(sqrt(optimValues.resnorm)*10000),par']);
        end
    end
    
    % Save updated result values:
    obj.pricerSettings.optim = res_obj;
    
end
function p = ModifyAndGetPrices(obj,p_fun,parNames,par,parIdx,parEqual,fSet,fInParNames,...
                                fOutPropNames)
%
%    Auxiliary function used by the 'Calibrate' method. Modifies model and thereafter 
%    computes prices. 
%
% -----------------------------------------------------------------------------------------------
%  Parameters  
% -----------------------------------------------------------------------------------------------
%   p_fun:           [1x1 function]   A function handle with no arguments returning the prices 
%                                     of some contracts under the model. Function call should 
%                                     refer to the object itself.
%   parNames:        [1xN cell]       See the 'Calibrate' method.
%   par:             [Mx1 real]       New parameter values.
%   parIdx:          [Mx1 integer]    Indices relating each element of the 'par' vector to a  
%                                     an element of 'parNames'.
%   parEqual:        [1xK cell]       See the 'Calibrate' method.
%   fSet:            [1xL cell]       See the 'Calibrate' method.
%   fInParNames:     [1xL cell]       See the 'Calibrate' method.
%   fOutPropNames:   [1xL cell]       See the 'Calibrate' method.
%
% -----------------------------------------------------------------------------------------------
%  Output
% -----------------------------------------------------------------------------------------------
%   p:               [Np x 1 real]    Prices or implied volatility. 
%
% -----------------------------------------------------------------------------------------------
    
    % Update parameters in object:
    for i=1:size(parNames,2)
        if isprop(obj,parNames{i})
            if isnumeric(obj.(parNames{i})) && max(size(isnumeric(obj.(parNames{i})))-[1,1])==0
                obj.(parNames{i}) = par(parIdx==i);
            elseif strcmpi(class(obj.(parNames{i})),'CurveClass')
                obj.(parNames{i}).values = par(parIdx==i);
            end
        end
    end
    
    % Update additional fixed parameters:
    for i=1:size(parEqual,2)
        obj.(parEqual{i}{1}) = obj.(parEqual{i}{2});
    end
    
    % Additional function calls:
    for i=1:size(fSet,2)
        fTmp = fSet{i};
        inNames = fInParNames{i};
        outNames = fOutPropNames{i};
        extraInputs = cell(size(inNames));
        for j=1:size(inNames,2)
            idxTmp = find(ismember(parNames,inNames{j})) == parIdx;
            extraInputs{j} = par(idxTmp);
        end
        outputTmp = fTmp(obj,extraInputs{:});
        for j=1:size(outNames,2)
            if isscalar(outputTmp)
                obj.(outNames{j}) = outputTmp;
            else
                obj.(outNames{j}) = outputTmp.DeepCopy();
            end
        end
    end

    % Compute prices:
    p = p_fun();
    
end
end
end

