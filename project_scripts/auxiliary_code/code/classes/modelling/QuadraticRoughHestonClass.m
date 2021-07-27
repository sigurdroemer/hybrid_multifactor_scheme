classdef QuadraticRoughHestonClass < PricingModelClass & handle
%
%   Implements the quadratic rough Heston model of (Gatheral et al., 2020). Let r(t) and 
%   delta(t) respectively denote the risk-free interest rate and continuous dividends yield. 
%   Both are assumed deterministic functions of time. Let (W_1(t),W_perp(t)) be a set of 
%   independent Brownian motions and define W_2(t) = rho*W_1(t) + sqrt(1-rho^2)*W_perp(t), 
%   rho in [-1,1]. The asset price S(t) is then assumed to have risk-neutral dynamics
%
%       dS(t) = S(t)(r(t) - delta(t))dt + S(t)sqrt(V(t))dW_2(t)
%
%   where 
%
%      V(t) = a*(Z(t)-b)^2 + c
%
%      Z(t) = zeta(t) + int_0^t K(t-s) eta*sqrt(V(s)) dW_1(s) 
%
%   with zeta and K being deterministic functions (K completely monotone).
% 
% ------------------------------------------------------------------------------------------------
%  Properties  
% ------------------------------------------------------------------------------------------------
%   a:        [1x1 real]                  Sensitivity of volatility to feedback of returns.
%   b:        [1x1 real]                  Assymmetry of feedback effect.
%   c:        [1x1 real]                  Minimum instantaneous variance.
%   eta:      [1x1 real]                  Volatility of Z(t).
%   rho:      [1x1 real]                  Correlation parameter.
%   K:        [1x1 KernelFunctionClass]   Kernel function. Must be vectorised.
%   zeta:     [1x1 CurveClass]            Mean value term structure of Z(t) process.
%
%   More properties are inherited from the PricingModelClass. The most important being the 
%   obj.pricerSettings property where the settings for the pricing algorithm are set. Since only 
%   Monte Carlo pricing is implemented for this class this property must be of the 
%   MonteCarloPricerSettingsClass type. You should consult the description of that class for more 
%   details on the possible settings. A few additional restrictons apply on top of those explained
%   in that class:
%
%       o obj.pricerSettings.simulation.scheme: Only options are 'hybrid_multifactor'.
%       o obj.pricerSettings.simulation.kernel_approximation.kappa: Can be any non-negative 
%         integer except if VIX2 needs to be simulated (e.g. if VIX options are requested). 
%         Then it must be set to 1.
%       o obj.pricerSettings.price_estimation.S.control_variate: Supported values are 'none' and 
%         'asset_price'. The latter uses S(T) as a control variate if conditional Monte Carlo
%         is not used. If conditional Monte Carlo is used, then we use S1(T) as a control variate.
%         See the 'simulate' method for a definition of the stochastic process S1(T).
%       o obj.pricerSettings.price_estimation.VIX.conditional_monte_carlo: Must be set to false.
%       o obj.pricerSettings.price_estimation.VIX.control_variate: Supported values are 'none' and 
%         'VIX2'. The latter uses VIX(T)^2 as a control variate.
%
% ------------------------------------------------------------------------------------------------
%  References
% ------------------------------------------------------------------------------------------------
%   o Jim Gatheral, Paul Jusselin and Mathieu Rosenbaum. The Quadratic Rough Heston Model and 
%     the Joint S&P 500/VIX Smile Calibration Problem, 2020. Working paper available at 
%     ssrn.com/abstract_id=3514894.
%
% ------------------------------------------------------------------------------------------------

properties        
    a
    b
    c    
    eta
    rho
    K
    zeta
end
    
methods
   function obj = QuadraticRoughHestonClass(varargin)
    % 	
    %   Constructor.
    %
    
        % Set object properties:
        obj.ParseConstructorInputs(varargin{:});
        
        % Set default settings for pricing algorithm:
        if isempty(obj.pricerSettings)
            % Price estimation settings:
            priceEst = struct;
            priceEst.antithetic = false;
            
            priceEst.S.control_variate = 'asset_price';
            priceEst.S.conditional_monte_carlo = false;
            priceEst.S.option_type = 'otm';
            priceEst.VIX.control_variate = 'none';
            priceEst.VIX.conditional_monte_carlo = false;
            priceEst.VIX.option_type = 'otm';

            % Simulation settings:
            sim = struct;
            sim.scheme = 'hybrid_multifactor';
            sim.explicit = false;
            
            KSettings = struct;
            KSettings.kappa = 1;
            KSettings.n_multiplier = 1;
            KSettings.n_cap = 15000;
            KSettings.error_measure = 'uniform_relative';
            KSettings.epsilon = 10^(-3);
            sim.kernel_approximation = KSettings;
            
            n = [50000;25000;12500;6400;3200;500];
            tn = [0.004;0.008;0.016;0.032;0.2;Inf];
            n_vix = [256;128;128;128;64;64];

            obj.pricerSettings = MonteCarloPricerSettingsClass('n',n,'tn',tn,...
                                                               'n_VIX2_integrate',n_vix,...
                                                               'n_EVIX2_integrate',n_vix,...
                                                               'N',50000,...
                                                               'price_estimation',priceEst,...
                                                               'simulation',sim,'numRepeat',1);
        else
            if isempty(obj.pricerSettings.tn) && size(obj.pricerSettings.n,1) == 1
                obj.pricerSettings.tn = Inf;
            end
        end

        % Set or adjust the zeta curve if needed:
        if isprop(obj,'zeta') && isempty(obj.zeta)
            obj.zeta = CurveClass('gridpoints',1,'values',0);
        end
        if isprop(obj,'zeta') && isnumeric(obj.zeta)
            obj.zeta = CurveClass('gridpoints',1,'values',obj.zeta);
        end
        
        % Adjust the yield and dividend yield curve if needed:
        if isprop(obj,'y') && isnumeric(obj.y)
            obj.y = CurveClass('gridpoints',1,'values',obj.y);
        end        
        if isprop(obj,'q') && isnumeric(obj.q)
            obj.q = CurveClass('gridpoints',1,'values',obj.q);
        end        
        
   end
   function [pS,idxCallS,seS,KS_out,kS_out,pV,idxCallV,seV,KV_out,kV_out,futV,seFutV,ivIdx] ...
                              = GetPricesSub(obj,KS,kS,qS,ttmS,cartProdS,KV,kV,qV,ttmV,cartProdV)
   %
   %    Computes prices of put and call options. The function is not intended to be used by an 
   %    end-user but only as an auxiliary function to be used by the 'GetPrices' method of the 
   %    PricingModelClass.
   %
   %    Remark: We assume that all input expirations use the same number of steps per year 
   %    (according to the obj.pricerSettings property).
   %
   % ---------------------------------------------------------------------------------------------
   %  Parameters
   % ---------------------------------------------------------------------------------------------
   % *  Exactly one of the inputs KS, kS and qS must be specified; the other should be left empty.
   %    The same goes for the inputs KV, kV and qV.
   %
   %    KS:             [Ns x 1 real] Strikes for options on S.   
   %    kS:             [Ns x 1 real] Log-moneyness for options on S.
   %    qS:             [Ns x 1 real] Quantiles for strikes on S.
   %    ttmS:           [Ms x 1 real] Expirations for options on S.
   %    cartProdS:      [1x1 logical] If true then we return prices for the cartesian product of  
   %                                  the moneyness and maturity vectors for options on S. Else we 
   %                                  assume Ns = Ms and return prices row by row.
   %    KV:             [Nv x 1 real] Strikes for options on VIX(t).
   %    kV:             [Nv x 1 real] Log-moneyness for options on VIX(t).
   %    qV:             [Nv x 1 real] Quantiles for strikes on VIX(t).
   %    ttmV:           [Mv x 1 real] Expirations for options VIX(t).
   %    cartProdV:      [1x1 logical] As for 'cartProdS' but for the VIX options.
   %  
   %  If, say, options on S(t) are not desired then the corresponding parameters can entirely be 
   %  left empty. Same goes if options on VIX(t) are not desired. 
   %
   % ---------------------------------------------------------------------------------------------    
   %  Outputs
   % ---------------------------------------------------------------------------------------------
   %    pS:         [Ns x 1 or Ns x Ms real]    Prices of either calls or puts on S.
   %    idxCallS:   [Ns x 1 or Ns x Ms logical] True ~ call, false ~ put for options on S.
   %    seS:        [Ns x 1 or Ns x Ms real]    Standard errors for options on S.
   %    KS_out:     [Ns x 1 or Ns x Ms real]    Strikes for options on S.   
   %    kS_out:     [Ns x 1 or Ns x Ms real]    Log-moneyness for options on S.
   %    pV:         [Nv x 1 or Nv x Mv real]    Prices of either call or put VIX options.
   %    idxCallV:   [Nv x 1 or Nv x Mv logical] True ~ call, false ~ put for VIX options.
   %    seV:        [Nv x 1 or Nv x Mv real]    Standard errors for VIX options.
   %    KV_out:     [Nv x 1 or Nv x Mv real]    Strikes for VIX options.
   %    kV_out:     [Nv x 1 or Nv x Mv real]    Log-moneyness for VIX options.
   %    futV:       [Nv x 1 or Nv x Mv real]    VIX futures prices.
   %    seFutV:     [Nv x 1 or Nv x Mv real]    Standard errors for VIX futures prices.
   %    ivIdx:      [1x1 logical]               True ~ prices are returned as implied volatility.
   %
   % ---------------------------------------------------------------------------------------------   
        
       % Initialize outputs:
       [pS,idxCallS,seS,KS_out,kS_out,pV,idxCallV,seV,KV_out,kV_out,futV,seFutV,ivIdx] = deal([]);
       
       % Note down what type of options are requested:
       doPriceS = ~isempty(ttmS);
       doPriceVIX = ~isempty(ttmV);
       
       if ~doPriceS && ~doPriceVIX
           return;
       end
       
       % Keep track of how moneyness is inputted:
       S_uses_q = ~isempty(qS);
       S_uses_k = ~isempty(kS);
       S_uses_K = ~isempty(KS);
       Ns = max([size(qS,1),size(kS,1),size(KS,1)]);
       VIX_uses_q = ~isempty(qV);
       VIX_uses_k = ~isempty(kV);
       VIX_uses_K = ~isempty(KV);
       Nv = max([size(qV,1),size(kV,1),size(KV,1)]);
       
       % Extract pricer settings:
       settings = obj.pricerSettings;
       optTypeS = settings.price_estimation.S.option_type;
       optTypeV = settings.price_estimation.VIX.option_type;
       cvS = settings.price_estimation.S.control_variate;
       cvV = settings.price_estimation.VIX.control_variate;
       condMC_S = settings.price_estimation.S.conditional_monte_carlo;
       condMC_V = settings.price_estimation.VIX.conditional_monte_carlo;
       antithetic = settings.price_estimation.antithetic;
       prec = settings.precision;
       
       if condMC_V
           error(['QuadraticRoughHestonClass:GetPricesSub: Conditional Monte Carlo is not ',...
                  'supported for VIX options.']);
       end
       if ~(strcmpi(cvS,'none') || strcmpi(cvS,'asset_price'))
           error(['QuadraticRoughHestonClass:GetPricesSub: ',...
                  'obj.price_estimation.S.control_variate ',...
                  'only supports options ''none'' and ''asset_price''.']);
       end
       if ~(strcmpi(cvV,'none') || strcmpi(cvV,'VIX2'))
           error(['QuadraticRoughHestonClass:GetPricesSub: ',...
                  'obj.price_estimation.VIX.control_variate ',...
                  'only supports options ''none'' and ''VIX2''.']);
       end     
       
       % Determine number of steps for simulation:
       [idxGrp,nStepsActual] = settings.GetNumSteps(unique([ttmS;ttmV]));
       
       % Determine number of integration points for VIX integrations:
       nI_VIX2 = settings.n_VIX2_integrate(idxGrp(1));
       nI_EVIX2 = settings.n_EVIX2_integrate(idxGrp(1));

       % Set extra parameters/settings for the simulation part:
       extraInput = {'rndNumbers',settings.random,'antithetic',antithetic,'n_is','total',...
                     'nI_VIX2',nI_VIX2};

       % Determine simulation variables and time points:
       ttmS_adj = sort(unique(ttmS))';
       ttmV_adj = sort(unique(ttmV))';
       [sim_variables, tOut] = deal({});
       if doPriceS
           if ~condMC_S
              sim_variables = {sim_variables{:},'S'};
              tOut = {tOut{:},ttmS_adj};
           else
              sim_variables = {sim_variables{:},'S1','QV'};
              tOut = {tOut{:},ttmS_adj,ttmS_adj};
              if ~isempty(qS)
                  sim_variables = {sim_variables{:},'S'};
                  tOut = {tOut{:},ttmS_adj};
              end
           end
           
       end
       if doPriceVIX
           sim_variables = {sim_variables{:},'VIX2'};
           tOut = {tOut{:},ttmV_adj};
       end
       
       % Simulate paths:
       paths0 = obj.Simulate(settings.N,max(nStepsActual),sim_variables,tOut,extraInput{:});
       
       % Price options on S:
       if doPriceS
            % Initialise outputs:
            tt_S = tOut{1};
            if cartProdS
                [pS,idxCallS,seS,kS_out,KS_out] = deal(zeros(Ns,size(tt_S,2)));                
            else
                [pS,idxCallS,seS,kS_out,KS_out] = deal(zeros(Ns,1));
            end

            % Compute prices for each expiration separately:
            for i=1:size(tOut{1},2)
                % First adjust 'paths0' to fit requirements of 'MonteCarloEstimation' method:
                pathsS = struct;
                pathsS.t = tt_S(i);
                if condMC_S
                    pathsS.S1 = paths0.S1.values(:,i);                       
                    pathsS.QV = paths0.QV.values(:,i);      
                    if S_uses_q
                        pathsS.S = paths0.S.values(:,i);
                    end
                else
                    pathsS.S = paths0.S.values(:,i);                       
                end

                % Compute strikes from quantiles:
                F = obj.s0.*exp((obj.y.Eval(tt_S(i))-obj.q.Eval(tt_S(i))).*tt_S(i));
                if cartProdS
                    if S_uses_q
                        KS_out(:,i) = F*quantile(pathsS.S,qS.',1);
                        kS_out(:,i) = log(KS_out(:,i)./F);
                    elseif S_uses_k
                        kS_out(:,i) = kS;
                        KS_out(:,i) = F.*exp(kS);
                    elseif S_uses_K
                        KS_out(:,i) = KS;
                        kS_out(:,i) = log(KS./F);
                    end
                else
                    idxT = ttmS == tt_S(i);
                    if S_uses_q
                        KS_out(idxT) = F*quantile(pathsS.S,qS(idxT).',1).';
                        kS_out(idxT) = log(KS_out(idxT)./F);
                    elseif S_uses_k
                        kS_out(idxT) = kS(idxT);
                        KS_out(idxT) = F.*exp(kS(idxT));
                    elseif S_uses_K
                        KS_out(idxT) = KS(idxT);
                        kS_out(idxT) = log(KS(idxT)./F);
                    end
                end

                % Estimate option prices:
                if cartProdS
                    [pS(:,i),idxCallS(:,i),seS(:,i)] = obj.MonteCarloEstimation(kS_out(:,i),...
                                                  tt_S(i),false,pathsS,...
                                                  optTypeS,condMC_S,cvS,antithetic,prec,F,...
                                                  obj.y.DeepCopy(),obj.rho);
                else
                    [pS(idxT),idxCallS(idxT),seS(idxT)] = obj.MonteCarloEstimation(...
                                                  kS_out(idxT),...
                                                  tt_S(i),false,pathsS,...
                                                  optTypeS,condMC_S,cvS,antithetic,prec,F,...
                                                  obj.y.DeepCopy(),obj.rho);                        
                end

            end            
       end
       
       % Price VIX options:
       if doPriceVIX
            % Initialise outputs:
            tt_VIX = tOut{end};
            if cartProdV
                [pV,idxCallV,seV,kV_out,KV_out,futV,seFutV] = deal(zeros(Nv,size(tt_VIX,2))); 
            else
                [pV,idxCallV,seV,kV_out,KV_out,futV,seFutV] = deal(zeros(Nv,1));
            end
           
            % Compute prices for each expiration separately:
            for i=1:size(tt_VIX,2)
               % Reset path objects:
               pathsV = [];

               % First adjust 'paths0' to fit requirements of 'MonteCarloEstimation' method:
               pathsV.t = tt_VIX(i);
               pathsV.S = sqrt(paths0.VIX2.values(:,i));

               % Add control variate:
               if strcmpi(cvV,'VIX2')
                    pathsV.Y = paths0.VIX2.values(:,i);
                    pathsV.EY = obj.GetExpectedValueVIX2(tt_VIX(i),nI_EVIX2);
                    cvMCEst = 'explicit';
               else
                    cvMCEst = cvV;
               end

               % Compute VIX futures:
               futV_tmp = mean(pathsV.S);
               se_futV_tmp = std(pathsV.S)./sqrt(size(pathsV.S,1));
               if cartProdV
                    futV(:,i) = futV_tmp;
                    seFutV(:,i) = se_futV_tmp;
               else
                   idxT = ttmV == tt_VIX(i);
                   futV(idxT) = futV_tmp;
                   seFutV(idxT) = se_futV_tmp;
               end
               pathsV.S = pathsV.S ./ futV_tmp;

               if cartProdV
                   if VIX_uses_q
                        % Compute strikes from quantiles:
                        KV_out(:,i) = futV_tmp*quantile(pathsV.S,qV.',1);
                        kV_out(:,i) = log(KV_out(:,i)./futV_tmp);
                   elseif VIX_uses_k
                        kV_out(:,i) = kV;
                        KV_out(:,i) = futV_tmp*exp(kV);
                   elseif VIX_uses_K
                        KV_out(:,i) = KV;
                        kV_out(:,i) = log(KV./futV_tmp);
                   end
               else
                   if VIX_uses_q
                        % Compute strikes from quantiles:
                        KV_out(idxT) = futV_tmp*quantile(pathsV.S,qV(idxT).',1).';
                        kV_out(idxT) = log(KV_out(idxT)./futV_tmp);
                   elseif VIX_uses_k
                        kV_out(idxT) = kV(idxT);
                        KV_out(idxT) = futV_tmp*exp(kV(idxT));
                   elseif VIX_uses_K
                        KV_out(idxT) = KV(idxT);
                        kV_out(idxT) = log(KV(idxT)./futV_tmp);
                   end
               end

               % Being careful we here reuse the 'MonteCarloEstimation' method which strictly
               % speaking is only implemented for options on S:
               if cartProdV
                   [pV(:,i),idxCallV(:,i),seV(:,i)] = obj.MonteCarloEstimation(...
                                                         kV_out(:,i),tt_VIX(i),...
                                                         cartProdV,pathsV,optTypeV,...
                                                         condMC_V,cvMCEst,antithetic,...
                                                         prec,futV_tmp,obj.y.DeepCopy(),[]);
               else
                   [pV(idxT),idxCallV(idxT),seV(idxT)] = obj.MonteCarloEstimation(...
                                                         kV_out(idxT),tt_VIX(i),...
                                                         cartProdV,pathsV,optTypeV,...
                                                         condMC_V,cvMCEst,antithetic,...
                                                         prec,futV_tmp,obj.y.DeepCopy(),[]);
               end
            end
            
       end

       % Prices are always returned in monetary terms:
       ivIdx = false;
       
   end
   function varargout = GenNumbersForSim(obj,N,M,anti,precision,onlyRetSpecs,variables)
   %
   %    Generates random numbers for method 'Simulation'.
   %
   % ---------------------------------------------------------------------------------------------
   %  Main parameters  
   % ---------------------------------------------------------------------------------------------
   %    N:   [1x1 integer]      Total number of paths.
   %    M:   [1x1 integer]      Total number of steps.
   %
   % ---------------------------------------------------------------------------------------------
   %  Optional parameters  
   % ---------------------------------------------------------------------------------------------
   %    anti:         [1x1 logical] Set to true to include antithetic samples (N must be divisible
   %                                by 2 in this case). Default is false.
   %    precision:    [1x1 string]  Options are 'single' and 'double' (default).
   %    onlyRetSpecs: [1x1 logical] If true we only return the variable names and the sizes of
   %                                each set of random numbers. See also the description under 
   %                                'Outputs'. Default value is false.   
   %    variables:    [1xL cell]    Names of variables that you wish to simulate. Must be a subset
   %                                of those available as described under method 'Simulate'.
   %
   % ---------------------------------------------------------------------------------------------
   %  Outputs
   % ---------------------------------------------------------------------------------------------
   %    If onlyRetSpecs = false the output is an object of class 'RandomNumbersClass' where the 
   %    'numbers' property is a struct with the following members:
   %
   %        o W1_Gaussians:     [N*M x (1+kappa) real]   i.i.d. N(0,1)'s related to W_1 Brownian 
   %                                                     motion. Output dimension depends on scheme.
   %        o Wperp_Gaussians:  [NxM real]               i.i.d. N(0,1)'s for W_perp Brownian motion.
   %
   %    If onlyRetSpecs = true there will be two ouputs:
   % 
   %        (1) [1x2 cell] The names of the members/variables that would be in the 
   %                       'RandomNumbersClass' if we had onlyRetSpecs = false. 
   %
   %        (2) [1x2 cell] Each element consists of a [1x2 integer] vector giving the dimension 
   %                       of each variable from (1) if we had onlyRetSpecs = false.
   %
   % ---------------------------------------------------------------------------------------------

       if ~exist('anti','var') || isempty(anti)
           anti = false;
       end
       if ~exist('onlyRetSpecs','var') || isempty(onlyRetSpecs)
           onlyRetSpecs = false;
       end
       if ~exist('precision','var') || isempty(precision)
           precision = 'double';
       end       
       
       scheme = obj.pricerSettings.simulation.scheme;
       if ~strcmpi(scheme,'hybrid_multifactor')
           error(['QuadraticRoughHestonClass:GetNumbersForSim: Scheme ', scheme, ...
               ' is not supported.']);
       end
       
       kappa = obj.pricerSettings.simulation.kernel_approximation.kappa;
       if kappa ~= 1 && any(ismember(variables,'VIX2'))
           error(['QuadraticRoughHestonClass:GetNumbersForSim: We currently require ',...
                  'kappa = 1 when simulating VIX2.']);
       end

       if anti
           % Check divisibility by 2
           if rem(N,2) ~= 0
               error(['QuadraticRoughHestonClass:GenNumbersForSim: When ',...
                      '''anti'' = true we must have N divisible by 2.']);
           end
       end
       
       if onlyRetSpecs
          if any(ismember({'S'},variables))
              vars = {'W1_Gaussians','Wperp_Gaussians'};
              sizes = {[N*M,kappa+1],[N,M]};
          else
              vars = {'W1_Gaussians'};  
              sizes = {[N*M,kappa+1]};
          end
          varargout = {vars,sizes};
          return;
       end

       % Generate numbers
       num = struct;
       if any(ismember({'S'},variables))
            num.Wperp_Gaussians = randn(N,M,precision);
       end

       if ~anti
            num.W1_Gaussians = randn(N*M,kappa+1,precision);
            
       else
           % Function 'HybridMultifactorScheme' interprets the input Gaussians such that the 
           % first N rows relate to step 1 for all paths, the next N rows to step 2 for all 
           % paths, etc. We therefore need to be extra careful when specifying the antithetic 
           % samples:

           idxOrig = false(N*M,1);
           iter = 0;
           for i=1:M
               for j=1:N
                   iter = iter + 1;
                   if j <= N/2
                      idxOrig(iter) = true;
                   end
               end
           end

           num.W1_Gaussians = zeros(N*M,kappa+1,precision);
           num.W1_Gaussians(idxOrig,:) = randn(N*M/2,kappa+1,precision);
           num.W1_Gaussians(~idxOrig,:) = -num.W1_Gaussians(idxOrig,:);     
           
       end

       % Create object:
       varargout = {RandomNumbersClass('numbers',num)};

   end   
   function GenNumbersForMC(obj,ttm,priceS,priceVIX,qS,seed)
   %
   %    Generates random numbers for one or several runs of the Monte Carlo pricer and stores 
   %    them in the obj.pricerSettings.random property.
   %
   % ---------------------------------------------------------------------------------------------
   %  Parameters  
   % ---------------------------------------------------------------------------------------------
   %    ttm:        [Nx1 real]      Expiries we wish to price using Monte Carlo.
   %    priceS:     [1x1 logical]   Set to true if you want to price option on S(t).
   %    priceVIX:   [1x1 logical]   Set to true if you want to price options on VIX(t).
   %    qS:         [1x1 logical]   Set to true if strikes for options on S(t) are specified
   %                                as quantiles of the terminal distribution.
   %    seed:       [1x1 integer    If non-empty we call rng(seed) before generating random numbers.
   %                 ... or empty]    
   % ---------------------------------------------------------------------------------------------

       settings = obj.pricerSettings;
       if ~isa(settings,'MonteCarloPricerSettingsClass')
           error(['QuadraticRoughHestonClass:GenNumbersForMC: Pricer settings must ', ...
                  'be an object of type ''MonteCarloPricerSettingsClass''.']);
       end

        % Find the total number of steps needed:
        [~,nTotalSteps] = settings.GetNumSteps(ttm);
        
        % Determine the required simulation variables:
        variables = {};
        if priceS
            if settings.price_estimation.S.conditional_monte_carlo
                variables = {variables{:},'S1','QV'};
                if qS
                    % We need S(t) to compute strikes from quantiles:
                    variables = {variables{:},'S'};
                end
            else
                variables = {variables{:},'S'};
            end
        end
        if priceVIX
            variables = {variables{:},'VIX2'};
        end
        
        % Generate random numbers:
        if ~isempty(seed);rng(seed);end
        settings.random = obj.GenNumbersForSim(settings.N,max(max(nTotalSteps)),...
                                               settings.price_estimation.antithetic,...
                                               settings.precision,false,variables);

   end
   function paths = Simulate(obj,N,n,variables,timepoints,varargin)
   %
   %    Simulates the model assuming an initial asset price of 1 and zero drift.
   %
   % ---------------------------------------------------------------------------------------------
   %  Main parameters  
   % ---------------------------------------------------------------------------------------------
   %   N:          [1x1 integer]   Number of paths to simulate.
   %   n:          [1x1 integer]   Number of time steps. By default is interpreted as number
   %                               of steps 'per unit of time'. Parameter can also specify the 
   %                               total number of steps by appropriately using the optional 
   %                               parameter 'n_is' (see further down). 
   %   variables:  [1xL cell]      Array with the names (strings) of the stochastic processes 
   %                               to be outputted. Options are: 
   %
   %                                o S:          Value of underlying asset.
   %                                o S1:         Value of the following stochastic process:
   %                                                S1(t) = exp(rho*int_0^t sqrt(V(s))dW_1(s)
   %                                                            - 0.5*rho^2*int_0^t V(s) ds)
   %                                o V:          The instantaneous variance process.
   %                                o QV:         Quadratic variation of log-returns, i.e.
   %                                                QV(t) = int_0^t V(s) ds.
   %                                o VIX2:       The VIX index squared.
   %                                o Z:          The Z-process; see the class description.
   %                                o U:          Let c = (c_1,...,c_m),gamm = (gamm_1,...,gamm_m) 
   %                                              for some m. We then approximate the kernel 
   %                                              function K as
   %
   %                                                K(t) approx(=) sum_{j=1}^m c_j e^(-gamm_j*t)
   %
   %                                              for t in [max(dt,kappa*dt),T] where 
   %
   %                                              T := max(T_basic,T_VIX2)
   %
   %                                              with 
   %
   %                                              T_basic = "maximum requested timepoint for 
   %                                                         variables S,S1,V,QV,U"
   %                                              T_VIX2 = "maximum requested timepoint for 
   %                                                        variable VIX2 plus 1/12" 
   %
   %                                              and where dt is the step size.
   %
   %                                              We then return the factors
   %
   %                                              U_j(t) := int_0^t exp(-gamm_j*(t-s))
   %                                                                eta*sqrt(V(s))dW_1(s)
   %                                              
   %                                              for j=1,...,m.
   %
   %   timepoints: [1xL cell]      Time points at which each stochastic process needs to be 
   %                               outputted. Each entry must be a row vector, must be sorted
   %                               in ascending order and values must be distinct. 
   %                               Remark: If the time points do not fit into the equidistant 
   %                               simulation grid, the value at such time points will be 
   %                               approximated by the value at the grid point just below. 
   %                               Also, zero as a 1x1 scalar is not allowed.
   %
   % ---------------------------------------------------------------------------------------------
   %  Additional parameters  
   % ---------------------------------------------------------------------------------------------
   %   * The following parameters should be given as name-value pairs. They are all optional
   %     except for 'nI_VIX2' which is a required parameter when variable 'VIX2' is requested.
   %
   %   nI_VIX2:     [1x1 integer]            Number of integration points for computing the VIX
   %                                         index (squared). We here use a trapezoidal rule. 
   %                                         Parameter simultaneously sets the number of points to 
   %                                         be used for solving the required Volterra equation to 
   %                                         obtain the forward variances.
   %
   %   n_is:       [1x1 string]              Specifies how to interpret the 'n' input. Options are
   %                                         'per_unit_of_time' and 'total' (i.e. total number 
   %                                         of time points till the largest timepoint in the
   %                                         'timepoints' input). Default is 'per_unit_of_time'.
   %
   %   antithetic: [1x1 logical]             If true we also return antithetic sample paths where 
   %                                         the values of the W_1(t) Brownian motion are flipped 
   %                                         to -W_1(t) (still N paths in total). It is then the 
   %                                         first N/2 paths returned that will be the original 
   %                                         paths and the bottom half that will be antithetic. 
   %                                         Default is antithetic = false. 
   %
   %   rndNumbers: [1x1 RandomNumbersClass]  Object containing the random numbers that will be 
   %                                         used for the simulation. If unspecified we simulate 
   %                                         them as needed.
   %
   %                                         Important: One should be careful to set this parameter
   %                                         correctly to avoid simulating paths with incorrect 
   %                                         properties. The recommended choice is therefore to 
   %                                         leave it unspecified. Carefully read the code if you
   %                                         wish to use it.
   %
   % ---------------------------------------------------------------------------------------------
   %  Output
   % ---------------------------------------------------------------------------------------------
   %   paths:  [1x1 struct] First and foremost contains a field for each element in 'variables'. 
   %                        In each such field then is a struct with the following fields: 
   %                            t:           [1xnT real] Time points.
   %                            t_truncated: [1xnT real] Time points truncated to simulation grid.
   %                            values:      [NxnT or NxmxnT or NxMxnT real] Simulated values. 
   %                                         Actual dimensions depend on the particular variable.
   %    
   %                        * nT = size(timepoints{i},2), i=1,...,L, for the i'th variable.
   %
   %                        If the simulation variable is 'U' we in addition include fields 'c' 
   %                        ([mx1 real])  'gamm' ([mx1 real]) for the coefficients explained 
   %                        under the 'variables' parameter.
   %
   % ---------------------------------------------------------------------------------------------

       % Parse inputs:
       inP = inputParser;
       addParameter(inP,'nI_VIX2',[]);
       addParameter(inP,'n_is','per_unit_of_time');
       addParameter(inP,'antithetic',false);
       addParameter(inP,'rndNumbers',[]);
       parse(inP,varargin{:});
       v2struct(inP.Results);
       
       % Check output variables are valid:
       expectedOutVars = {'S','S1','V','QV','VIX2','U','Z'};
       if sum(ismember(variables,expectedOutVars)) ~= size(variables,2)
           error(['QuadraticRoughHestonClass:Simulate: Some of the chosen output variables ',...
                  'are not supported.']);
       end  
       
       % Check for duplicates:
       if size(unique(variables),2) ~= size(variables,2)
           error(['QuadraticRoughHestonClass:Simulate: Values in ''variables'' parameter ',...
                  'must be unique.']);
       end
       
       if size(variables,1) ~= size(timepoints,1) || size(variables,2) ~= size(timepoints,2)
           error(['QuadraticRoughHestonClass:Simulate: The inputs ''variables'' and ',...
                  '''timepoints'' must have the same dimensions']);
       end
       
       % Compute maximum time point:
       T = 0; % ~ horison to simulate till
       T_approx_K  = 0; % ~ horizon to approximate K(t) till
       for i=1:size(timepoints,2)
           if size(timepoints{i},1) > 1
               error(['QuadraticRoughHestonClass:Simulate: The elements of ''timepoints'' ',...
                      'must all be row vectors.']);
           end
           if ~issorted(timepoints{i},'ascend')
               error(['QuadraticRoughHestonClass:Simulate: Each element of ''timepoints'' ',...
                      'must be sorted in ascending order.']);
           end
           if any(timepoints{i} < 0)
               error(['QuadraticRoughHestonClass:Simulate: The elements of ''timepoints'' ',...
                      'must be non-negative.']);
           end
           if size(timepoints{i},2) == 1 && timepoints{i} == 0
               error(['QuadraticRoughHestonClass:Simulate: The elements of ''timepoints'' ',...
                      'are not allowed to be scalar [1x1 real] with a zero value.']);
           end
           if size(unique(timepoints{i}),2) ~= size(timepoints{i},2)
               error(['QuadraticRoughHestonClass:Simulate: The elements of ''timepoints'' ',...
                      'must contain distinct values only.']);               
           end
           T = max([T,max(timepoints{i})]);
           if strcmpi(variables{i},'VIX2')
               T_approx_K = max([T_approx_K,max(timepoints{i})+1/12]);
           else
               T_approx_K = max([T_approx_K,max(timepoints{i})]);
           end
       end
       
       % Construct simulation grid: 
       if strcmpi(n_is,'per_unit_of_time')
            dt = 1/n;
            nStepsTotalSim = floor2(n*T);
       elseif strcmpi(n_is,'total')
            dt = T/n;
            nStepsTotalSim = n;
       else
            error(['QuadraticRoughHestonClass:Simulate: Invalid value of input ',...
                   'parameter ''n_is''.']);
       end
       n_per_unit = 1/dt;
       grid_pts = dt*(0:nStepsTotalSim);
       T = max(grid_pts);
       
       % Unpack a few settings:
       simSettings = obj.pricerSettings.simulation;
       prec = obj.pricerSettings.precision;
       
       % Extract model formulation:
       zeta0 = @(t)(obj.zeta.Eval(t));
       eta = obj.eta;rho = obj.rho;
       K = obj.K;
       sig = @(t,x)(eta*sqrt(obj.a.*(x-obj.b).^2+obj.c));
       
       % Note down which variables are requested:
       idxS = find(ismember(variables,'S'));
       idxS1 = find(ismember(variables,'S1'));
       idxV = find(ismember(variables,'V'));
       idxQV = find(ismember(variables,'QV'));
       idxU = find(ismember(variables,'U'));
       idxVIX2 = find(ismember(variables,'VIX2'));       
       idxZ = find(ismember(variables,'Z'));  
       
       if ~isempty(idxVIX2) && isempty(nI_VIX2)
           error(['QuadraticRoughHestonClass:Simulate: When simulating ''VIX2'' we require ',...
                  'the input parameter ''nI_VIX2'' to be set.']);
       end
       
       % Compute sum of exponentials approximation:
       kernelApproxSettings = obj.pricerSettings.simulation.kernel_approximation;
       kappa = kernelApproxSettings.kappa;
       if kappa ~= 1 && ~isempty(idxVIX2)
           error(['QuadraticRoughHestonClass:Simulate: We currently require ',...
                  'kappa = 1 when simulating VIX2.']);
       end    
       
       Npts = ceil2(n_per_unit*(T_approx_K-max(dt,kappa*dt))*kernelApproxSettings.n_multiplier);
       if Npts > kernelApproxSettings.n_cap
           Npts = kernelApproxSettings.n_cap;
           warning(['QuadraticRoughHestonClass:Simulate: Number of sampling points for ',...
                    'sum of exponentials approximation was truncated at ', ...
                    num2str(kernelApproxSettings.n_cap), ' points.']);
       end
       
       [c,gamm] = ApproximateKernel('BM2005','K',@(t)(obj.K.Eval(t)),'T',T_approx_K,...
                                    'delta',max(dt,dt*kappa),'N',Npts,...
                                    'epsilon',kernelApproxSettings.epsilon,...
                                    'error_measure',kernelApproxSettings.error_measure);

                                
       % Determine time points to initially return U-values for:
      if ~any(ismember(variables,{'U','VIX2'}))
          % No time points needed:
          tU_sim = [];
          returnU = false;
      else
          % We return all time points needed for the requested variables:
          idxVolRelatedVars = find(ismember(variables,{'U','VIX2'}));
          tU_sim = [];
          for i=1:size(idxVolRelatedVars,2)
              tU_sim = unique([tU_sim,timepoints{idxVolRelatedVars(i)}]);
          end
          returnU = true;
      end 
       
       % Determine time points to initially return V-values for:
       if ~isempty(idxS) || ~isempty(idxS1) || ~isempty(idxQV) 
           tV_sim = grid_pts;
       else
          idxVolVars = find(ismember(variables,{'V','U','VIX2','Z'}));
          tV_sim = [];
          for i=1:size(idxVolVars,2)
              tV_sim = unique([tV_sim,timepoints{idxVolVars(i)}]);
          end
       end

       % Simulate random variables:
       if isempty(rndNumbers)
            rndNumbers = obj.GenNumbersForSim(N,nStepsTotalSim,antithetic,prec,false,variables);
       end
       
       if any(ismember({'S','S1'},variables))
           returndW = true;
       else
           returndW = false;
       end
       
       % Simulate Z-process:
       switch simSettings.scheme
           case 'hybrid_multifactor'
               % Remark: Input nStepsTotalSim/T may not be integer valued as the description of 
               % 'HybridMultifactorScheme' requires it to be. However, a careful reading of
               % the code shows that non-integer values are also valid.
               [Zpaths,Upaths,dW1] = HybridMultifactorScheme(N,nStepsTotalSim/T,[],zeta0,0,sig,...
                                        gamm,c,kappa,'K',K,'returnU',returnU,...
                                        'positive',false,'tX',tV_sim,...
                                        'tU',tU_sim,'returndW',returndW,...
                                        'precision',prec,...
                                        'explicit',simSettings.explicit,...
                                        'Z',rndNumbers.numbers.W1_Gaussians(1:nStepsTotalSim*N,:));
               Upaths.c = c;
               Upaths.gamm = gamm;
               
           otherwise
               error('QuadraticRoughHestonClass:Simulate: Scheme is not supported.');
       end
       
       % Compute V process:
       Vpaths = struct;
       Vpaths.t = Zpaths.t;
       Vpaths.t_truncated = Zpaths.t_truncated;
       Vpaths.values = obj.a*(Zpaths.values-obj.b).^2 + obj.c;
       
       % Compute and save variables:
       paths = struct;
       
       % Z:
       if ~isempty(idxZ)
           idxt = LocateTimepointsInDiscretisation(Zpaths.t,timepoints{idxZ});
           paths.Z.t = timepoints{idxZ};
           paths.Z.t_truncated = Zpaths.t_truncated(idxt);
           paths.Z.values = Zpaths.values(:,idxt);
       end

       % V: 
       if ~isempty(idxV)
           idxt = LocateTimepointsInDiscretisation(Vpaths.t,timepoints{idxV});
           paths.V.t = timepoints{idxV};
           paths.V.t_truncated = Vpaths.t_truncated(idxt);
           paths.V.values = Vpaths.values(:,idxt);
       end
       
       % QV:
       if ~isempty(idxQV)
          QV = [zeros(size(Vpaths.values,1),1),cumsum(Vpaths.values(:,1:end-1),2)*dt];
          idxt = LocateTimepointsInDiscretisation(Vpaths.t,timepoints{idxQV});
          paths.QV.t = timepoints{idxQV};
          paths.QV.t_truncated = Vpaths.t_truncated(idxt);
          paths.QV.values = QV(:,idxt);
       end
       
       % U: 
       if ~isempty(idxU)
           idxt = LocateTimepointsInDiscretisation(Upaths.t,timepoints{idxU});
           paths.U.t = timepoints{idxU};
           paths.U.t_truncated = Upaths.t_truncated(idxt);
           paths.U.values = Upaths.values(:,:,idxt);      
           paths.U.c = Upaths.c;
           paths.U.gamm = Upaths.gamm;
       end
 
       % S1:
       if ~isempty(idxS1)
           dlogS1 = zeros(N,size(Vpaths.values,2),prec);
           dlogS1(:,2:end) = sqrt(Vpaths.values(:,1:end-1)).*rho.*dW1 ...
                            - 0.5*(rho^2)*Vpaths.values(:,1:end-1)*dt; 
           S1 = exp(cumsum(dlogS1,2));
           idxt = LocateTimepointsInDiscretisation(Vpaths.t,timepoints{idxS1});
           paths.S1.t = timepoints{idxS1};
           paths.S1.t_truncated = Vpaths.t_truncated(idxt);
           paths.S1.values = S1(:,idxt);
       end
              
       % S: 
       if ~isempty(idxS)
           dW2 = rho*dW1 + sqrt(1-rho^2)...
               *rndNumbers.numbers.Wperp_Gaussians(1:size(dW1,1),1:size(dW1,2))*sqrt(dt);               
           dlogS = zeros(N,size(Vpaths.values,2),prec);
           dlogS(:,2:end) = sqrt(Vpaths.values(:, 1:end-1)).*dW2(:,1:size(Vpaths.values,2)-1) ...
                            - 0.5 *Vpaths.values(:,1:end-1)*dt; 
           S = exp(cumsum(dlogS,2));
           idxt = LocateTimepointsInDiscretisation(Vpaths.t,timepoints{idxS});
           paths.S.t = timepoints{idxS};
           paths.S.t_truncated = Vpaths.t_truncated(idxt);
           paths.S.values = S(:,idxt);
       end
       
       % Remark: For the code below we currently require kappa = 1. That indeed kappa = 1, is 
       % verified further up in the code.
       
       % zeta:
       if ~isempty(idxVIX2)
           % Determine maturities needed for integration to get VIX2:
           DELTA = 1/12;
           dTau = DELTA/(nI_VIX2 - 1);
           tau = (0:dTau:DELTA);
           
           idxtU = LocateTimepointsInDiscretisation(Upaths.t,timepoints{idxVIX2});
           idxtV = LocateTimepointsInDiscretisation(Vpaths.t,timepoints{idxVIX2});
           
           paths.zeta.t = timepoints{idxVIX2};
           paths.zeta.t_truncated = Upaths.t_truncated(idxtU);
           paths.zeta.values = zeros(N,size(tau,2),size(idxtU,2));
           U = permute(Upaths.values(:,:,idxtU),[2,1,3]);
           for i=1:size(tau,2)
                if tau(i) >= dt
                    % Use sum of exponentials formulation:
                    if N == 1 || size(paths.zeta.t,2) == 1 
                        paths.zeta.values(:,i,:) = zeta0(paths.zeta.t + tau(i)) + squeeze(...
                                             sum(U.*(Upaths.c.*exp(-Upaths.gamm.*tau(i))),1)).';                        
                    else
                        paths.zeta.values(:,i,:) = zeta0(paths.zeta.t + tau(i)) + squeeze(...
                                             sum(U.*(Upaths.c.*exp(-Upaths.gamm.*tau(i))),1));
                    end
                else
                    % Use linear interpolation:
                    if N == 1 || size(paths.zeta.t,2) == 1
                        zeta_dt = zeta0(paths.zeta.t + dt) + squeeze(...
                                sum(U.*(Upaths.c.*exp(-Upaths.gamm.*dt)),1))';                        
                    else
                        zeta_dt = zeta0(paths.zeta.t + dt) + squeeze(...
                                sum(U.*(Upaths.c.*exp(-Upaths.gamm.*dt)),1));
                    end
                    frac = tau(i)/dt;
                    paths.zeta.values(:,i,:) = (1-frac)*Zpaths.values(:,idxtV) + frac*zeta_dt;
                end
           end
       end
       
       % VIX2:
       if ~isempty(idxVIX2)
            idxtU = LocateTimepointsInDiscretisation(Upaths.t,timepoints{idxVIX2});
            paths.VIX2.t = timepoints{idxVIX2};
            paths.VIX2.t_truncated = Upaths.t_truncated(idxtU);

            % Define squared kernel:
            [~,~,~,Ktilde] = ExpandProductOfKernels(obj.K,obj.K);
            
            % Approximate K-tilde:
            n_per_unit = 1/dTau;
            Npts = ceil2(n_per_unit*(DELTA-max(dt,kappa*dt))...
                        *kernelApproxSettings.n_multiplier);
            if Npts > kernelApproxSettings.n_cap
                Npts = kernelApproxSettings.n_cap;
                warning(['QuadraticRoughHestonClass:Simulate: Number of sampling ',...
                         'points for sum of exponentials approximation was truncated at ', ...
                         num2str(kernelApproxSettings.n_cap), ' points.']);
            end   
            [c_tilde,gamm_tilde] = ApproximateKernel('BM2005',...
                   'K',@(t)(Ktilde.Eval(t)),'T',DELTA,'delta',max(dTau,dTau*kappa),'N',Npts,...
                   'epsilon',kernelApproxSettings.epsilon,...
                   'error_measure',kernelApproxSettings.error_measure);     
            
            % Define f:
            ftu = obj.a*paths.zeta.values.^2 - 2*obj.a*obj.b*paths.zeta.values...
                  + obj.a*obj.b.^2 + obj.c;
               
            paths.VIX2.values = zeros(N,size(timepoints{idxVIX2},2));
            dummy = zeros(size(ftu,1),size(ftu,2)-1,kappa+1,prec);
            for i=1:size(paths.VIX2.t,2)
                % Run hybrid multifactor scheme:
                xi = HybridMultifactorScheme(N,n_per_unit,[],ftu(:,:,i),...
                                             @(t,x)((obj.a*obj.eta.^2)*x),0,...
                                             gamm_tilde,c_tilde,kappa,...
                                             'K',Ktilde,'positive',true,'tX',tau,...
                                             'explicit',simSettings.explicit,'W_bold',dummy);
     
                % Compute and save VIX2 values:
                paths.VIX2.values(:,i) = (1/DELTA)*(0.5*xi.values(:,1) ...
                                           + sum(xi.values(:,2:end-1),2) ...
                                           + 0.5*xi.values(:,end))*dTau;            

            end                      
       end
       
       % Clean 'paths' object before returning:
       if isfield(paths,'zeta')
           paths = rmfield(paths,'zeta');
       end

   end
   function zeta = GetZetaFromForwardVariance(obj,t,xi)
   %
   %    Computes E(Z(t)), t >= 0, for a given forward variance curve E(V(t)), t >= 0.
   %
   % ---------------------------------------------------------------------------------------------
   %  Parameters  
   % ---------------------------------------------------------------------------------------------       
   %    t:      [Nx1 real]      Time points.
   %    xi:     [1x1 function]  Function representing xi(t) := E(V(t)).
   %
   % ---------------------------------------------------------------------------------------------       
   %  Output
   % ---------------------------------------------------------------------------------------------       
   %    zeta:   [Nx1 real]  Values of E(Z(t)).
   %
   % ---------------------------------------------------------------------------------------------             
   
        K2_conv_xi = NaN(size(t));
        for i=1:size(t,1)
            K2_conv_xi(i) = obj.K.IntegrateVolterra(0,t(i),t(i),2,xi);
        end
        f = xi(t) - obj.a*(obj.eta^2)*K2_conv_xi;
        
        % Solve 2nd degree polynomial:
        a_ = obj.a;
        b_ = -2*obj.a*obj.b;
        c_ = obj.a*obj.b^2 + obj.c - f;
        zeta = (-b_ - sqrt(b_^2 - 4*a_*c_))./(2*a_);
   
   end
   function xi = GetForwardVariance(obj,t,dt,c_tilde,gamm_tilde)
   %
   %    Computes initial forward variances, that is, expectations of the form E(V_t|F_0).
   %
   % ---------------------------------------------------------------------------------------------
   %  Parameters  
   % ---------------------------------------------------------------------------------------------       
   %    t:            [Nx1 real]  Maturities.
   %    dt:           [1x1 real]  Step size for solving the Volterra equation.
   %    c_tilde:      [mx1 real]  Not intended for an end-user. Used by 'GetExpectedValueVIX2' for
   %                              performance reasons.
   %    gamm_tilde:   [mx1 real]  As for 'c_tilde'.
   %
   % ---------------------------------------------------------------------------------------------       
   %  Output
   % ---------------------------------------------------------------------------------------------       
   %    xi:     [Nx1 real] Forward variances.
   %
   % ---------------------------------------------------------------------------------------------       
      
       if size(t,2) > 1
           error('QuadraticRoughHestonClass:GetForwardVariance: t must be a column vector.');
       end
       
        % Extract model variables:
        a = obj.a;
        b = obj.b;
        eta = obj.eta;
        c = obj.c;
        
        % Define 'f' function:
        zeta = @(t)(obj.zeta.Eval(t));
        f0 = @(t)(a*zeta(t).^2 - 2*a*b*zeta(t) + a*b^2 + c);   

        % Create discretisation:
        n_per_unit = 1/dt;
        T = max(t);
        
        % Define squared kernel:
        [~,~,~,Ktilde] = ExpandProductOfKernels(obj.K,obj.K);
        
        % Approximate K-tilde:
        KSettings = obj.pricerSettings.simulation.kernel_approximation;
        kappa = KSettings.kappa;
        if (~exist('c_tilde','var') || isempty(c_tilde)) ...
                && (~exist('gamm_tilde','var') || isempty(gamm_tilde))
            Npts = ceil2(n_per_unit*(T-max(dt,dt*kappa))*KSettings.n_multiplier);
            if Npts > KSettings.n_cap
                Npts = KSettings.n_cap;
                warning(['QuadraticRoughHestonClass:GetForwardVariance: Number of sampling ',...
                         'points for sum of exponentials approximation was truncated at ', ...
                         num2str(KSettings.n_cap), ' points.']);
            end            
            [c_tilde,gamm_tilde] = ApproximateKernel('BM2005',...
                                                     'K',@(s)(Ktilde.Eval(s)),...
                                                     'T',T,'delta',max(dt,dt*kappa),'N',Npts,...
                                                     'epsilon',KSettings.epsilon,...
                                                     'error_measure',KSettings.error_measure);     
        else
            if size(c_tilde,2) > 1 || size(gamm_tilde,2) > 1 ...
                    || size(c_tilde,1) ~= size(gamm_tilde,1)
                
                error(['QuadraticRoughHestonClass:GetForwardVariance: ''c_tilde'' and ',...
                       '''gamm_tilde'' must be column vectors of the same size.']);
                
            end
        end
        
        % Solve equation by running hybrid multifactor scheme:
        dummy = zeros(1,floor2(max(t)*(1/dt)),kappa+1,obj.pricerSettings.precision);
        paths = HybridMultifactorScheme(1,1/dt,[],f0,@(t,x)((a*eta.^2)*x),0,...
                                        gamm_tilde,c_tilde,kappa,...
                                        'K',Ktilde,'positive',true,'tX',t',...
                                        'explicit',obj.pricerSettings.simulation.explicit,...
                                        'W_bold',dummy);
        xi = paths.values.';
        
   end
   function E_VIX2 = GetExpectedValueVIX2(obj,t,nI)
   %
   %    Computes expected values of VIX^2, i.e. expectations of the form E(VIX(t)^2|F_0).
   %
   %    Remarks: 
   %       o Function has not yet been optimised for large N.
   %
   % ---------------------------------------------------------------------------------------------
   %  Parameters  
   % ---------------------------------------------------------------------------------------------       
   %    t:      [Nx1 real]      Time points.
   %    nI:     [1x1 integer]   Number of integration points per month (i.e. per 1/12 unit 
   %                            of time).
   %
   % ---------------------------------------------------------------------------------------------       
   %  Output
   % ---------------------------------------------------------------------------------------------       
   %    E_VIX2: [Nx1 real] Expected values of VIX squared.
   %
   % ---------------------------------------------------------------------------------------------              
       
       if size(t,2) > 1
           error('QuadraticRoughHestonClass:GetExpectedValueVIX2: t must be a column vector.');
       end
   
       KSettings = obj.pricerSettings.simulation.kernel_approximation;
       kappa = KSettings.kappa;
       n_per_unit = 12*(nI-1); % nI - 1 =  number of steps per month
       dt = 1/n_per_unit; % step size will be adjusted for each element of t if needed
       DELTA = 1/12;
       
       % Compute distinct time points where we need forward variances:
       tt = [];
       ttCl = cell(size(t));
       for i=1:size(t,1)           
           tt0 = (t(i):dt:t(i)+DELTA)';
           if tt0(end) == t(i) + DELTA
               dtAdj = dt;
           else
               nSteps = size(tt0,1);
               dtAdj = DELTA/nSteps;
           end
           ttCl{i} = (t(i):dtAdj:t(i)+DELTA)';
           tt = [tt;ttCl{i}];
       end
       uniqT = unique(tt);
       T = max(uniqT);
       
       % Compute number of points for approximating kernel:
       Npts = ceil2(n_per_unit*(T-max(dt,dt*kappa))*KSettings.n_multiplier);
       if Npts > KSettings.n_cap
           Npts = KSettings.n_cap;
           warning(['QuadraticRoughHestonClass:GetExpectedValueVIX2: Number of sampling ',...
                    'points for sum of exponentials approximation was truncated at ', ...
                    num2str(KSettings.n_cap), ' points.']);
       end    

       % Approximate kernel:
       [c_tilde,gamm_tilde] = ApproximateKernel('BM2005',...
                                                'K',@(s)(obj.K.Eval(s).^2),...
                                                'T',T,'delta',max(dt,dt*kappa),'N',Npts,...
                                                'epsilon',KSettings.epsilon,...
                                                'error_measure',KSettings.error_measure);          

       % Compute forward variances:
       xi = obj.GetForwardVariance(uniqT,dt,c_tilde,gamm_tilde);                                        
                                        
       % Loop and compute E(VIX2) for each input:
       E_VIX2 = zeros(size(t));
       for i=1:size(t,1)
           % Find relevant values:
           idx = ismembertol(uniqT,ttCl{i});
           xiSub = xi(idx);
           
           % We use a trapezoidal rule:
           E_VIX2(i) = (1/DELTA)*(0.5*xiSub(1) + sum(xiSub(2:end-1)) + 0.5*xiSub(end))*dt;
       end
       
   end   
end

end
