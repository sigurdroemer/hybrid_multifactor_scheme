classdef MixedrBergomiClass < PricingModelClass & handle
%
%   Implements the following stochastic volatility model: Let r(t) and delta(t) respectively 
%   denote the risk-free interest rate and continuous dividends yield. Both are assumed 
%   deterministic functions of time. The asset price S(t) is assumed to have risk-neutral dynamics
%
%     dS(t) = S(t)(r(t) - delta(t))dt + S(t)sqrt(V(t))dW1(t)
%
%   where 
%
%     V(t) = c + zeta(t)*(theta(t)*exp(X_1(t)) + (1-theta(t))*exp(X_2(t)))
%
%   with
%
%   X_1(t) = eta_1*int_0^t K_1(t-s)dW2(t) - 0.5*eta_1^2*int_0^t K_1(s)^2 ds
%
%   X_2(t) = eta_2*int_0^t K_2(t-s)dW3(t) - 0.5*eta_2^2*int_0^t K_2(s)^2 ds
%
%   with zeta, theta, K_1, K_2 being deterministic functions (K's completely monotone) and 
%   where c, eta_1, eta_2 > 0, 0 <= theta(t) <= 1 and W = (W1,W2,W3) have a general correlation
%   structure as specified by parameters rho_12, rho_13, rho_23. 
% 
% -----------------------------------------------------------------------------------------------
%  Properties  
% -----------------------------------------------------------------------------------------------
%   c:        [1x1 real]                Minimum instantaneous variance.
%   eta_1:    [1x1 real]                Volatility of X_1.
%   eta_2:    [1x1 real]                Volatility of X_2.
%   rho_12:   [1x1 real]                Correlation between W1 and W2.
%   rho_13:   [1x1 real]                Correlation between W1 and W3.
%   rho_23:   [1x1 real]                Correlation between W2 and W3.
%   K_1:      [1x1 KernelFunctionClass] Kernel for X_1.
%   K_2:      [1x1 KernelFunctionClass] Kernel for X_2.
%   theta:    [1x1 Curveclass]          Mixing weight function.
%   zeta:     [1x1 CurveClass]          Curve to control the mean term structure of V(t). 
%                                       The forward variance curve xi(t) is related as 
%                                       xi(t) := zeta(t) + c.
%
%   More properties are inherited from the PricingModelClass. The most important being the 
%   obj.pricerSettings property where the settings for the pricing algorithm are set. Since only 
%   Monte Carlo pricing is implemented for this class this property must be of the 
%   MonteCarloPricerSettingsClass type. You should consult the description of that class for more 
%   details on the possible settings. A few additional restrictons apply on top of those explained
%   in that class:
%
%       o obj.pricerSettings.simulation.scheme: Only option is 'hybrid_multifactor'.
%       o obj.pricerSettings.simulation.kernel_approximation.kappa: Can be any non-negative 
%         integer except if VIX2 needs to be simulated (e.g. if VIX options are requested). 
%         Then it must be set to 1.
%       o obj.pricerSettings.price_estimation.S.control_variate: Supported values are 'none' and 
%         'asset_price'. The latter uses S(T) as a control variate.
%       o obj.pricerSettings.price_estimation.S.conditional_monte_carlo: Must be set to false.
%       o obj.pricerSettings.price_estimation.VIX.conditional_monte_carlo: Must be set to false.
%       o obj.pricerSettings.price_estimation.VIX.control_variate: Supported values are 'none' and 
%         'VIX2'. The latter uses VIX(T)^2 as a control variate.
%
% -----------------------------------------------------------------------------------------------

properties      
    c
    eta_1
    eta_2
    rho_12
    rho_13
    rho_23
    K_1
    K_2
    theta
    zeta
end
    
methods
   function obj = MixedrBergomiClass(varargin)
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
            n_vix = [32;32;32;32;32;32];

            obj.pricerSettings = MonteCarloPricerSettingsClass('n',n,'tn',tn,...
                                                               'n_VIX2_integrate',n_vix,...
                                                               'n_EVIX2_integrate',[],...
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
        
        % Set or adjust the theta curve if needed:
        if isprop(obj,'theta') && isempty(obj.theta)
            obj.theta = CurveClass('gridpoints',1,'values',0);
        end
        if isprop(obj,'theta') && isnumeric(obj.theta)
            obj.theta = CurveClass('gridpoints',1,'values',obj.theta);
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
       
       if condMC_S
           error(['MixedrBergomiClass:GetPricesSub: Conditional Monte Carlo is not ',...
                  'supported for options on S.']);
       end
       if condMC_V
           error(['MixedrBergomiClass:GetPricesSub: Conditional Monte Carlo is not ',...
                  'supported for VIX options.']);
       end
       if ~(strcmpi(cvS,'none') || strcmpi(cvS,'asset_price'))
           error(['MixedrBergomiClass:GetPricesSub: obj.price_estimation.S.control_variate ',...
                  'only supports options ''none'' and ''asset_price''.']);
       end
       if ~(strcmpi(cvV,'none') || strcmpi(cvV,'VIX2'))
           error(['MixedrBergomiClass:GetPricesSub: obj.price_estimation.VIX.control_variate ',...
                  'only supports options ''none'' and ''VIX2''.']);
       end

       % Determine number of steps for simulation:
       [idxGrp,nStepsActual] = settings.GetNumSteps(unique([ttmS;ttmV]));
       
       % Determine number of integration points for VIX integrations:
       nI_VIX2 = settings.n_VIX2_integrate(idxGrp(1));

       % Set extra parameters/settings for the simulation part:
       extraInput = {'rndNumbers',settings.random,'antithetic',antithetic,'n_is','total',...
                     'nI_VIX2',nI_VIX2};

       % Determine simulation variables and time points:
       ttmS_adj = sort(unique(ttmS))';
       ttmV_adj = sort(unique(ttmV))';
       [sim_variables,tOut] = deal({});
       if doPriceS
          sim_variables = {sim_variables{:},'S'};
          tOut = {tOut{:},ttmS_adj};
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
                pathsS.S = paths0.S.values(:,i);                       

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
                                                  obj.y.DeepCopy(),[]);
                else
                    [pS(idxT),idxCallS(idxT),seS(idxT)] = obj.MonteCarloEstimation(...
                                                  kS_out(idxT),...
                                                  tt_S(i),false,pathsS,...
                                                  optTypeS,condMC_S,cvS,antithetic,prec,F,...
                                                  obj.y.DeepCopy(),[]);                        
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
                    pathsV.Y = paths0.VIX2.values;
                    pathsV.EY = obj.zeta.Integrate(paths0.VIX2.t',paths0.VIX2.t'+1/12)'*12 + obj.c;
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
   function varargout = GenNumbersForSim(obj,N,M,anti,precision,onlyRetSpecs)
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
   %
   % ---------------------------------------------------------------------------------------------
   %  Outputs
   % ---------------------------------------------------------------------------------------------
   %    If onlyRetSpecs = false the output is an object of class 'RandomNumbersClass' where the 
   %    'numbers' property is a struct with the following members:
   %
   %        o Z1_Gaussians:     [N*M x (1+2*kappa) real]   i.i.d. N(0,1)'s.
   %        o Z2_Gaussians:     [N*M x (1+2*kappa) real]   i.i.d. N(0,1)'s.
   %        o Z3_Gaussians:     [N*M x (1+kappa) real]     i.i.d. N(0,1)'s.
   %
   %    If onlyRetSpecs = true there will be two ouputs:
   % 
   %        (1) [1x2 cell] The names of the members/variables that would be in the 
   %                       'RandomNumbersClass' if we had onlyRetSpecs = false. 
   %
   %        (2) [1x2 cell] Each element consists of a [1x2 integer] vector giving the dimension 
   %                       of each variable from (1) if we had onlyRetSpecs = false.
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
           error(['MixedrBergomiClass:GetNumbersForSim: Scheme ', scheme, ' is not supported.']);
       end
       
       kappa = obj.pricerSettings.simulation.kernel_approximation.kappa;

       if anti
           % Check divisibility by 2
           if rem(N,2) ~= 0
               error(['MixedrBergomiClass:GenNumbersForSim: When ',...
                      '''anti'' = true we must have N divisible by 2.']);
           end
           nIndep = N/2; 
       else
           nIndep = N;
       end

       if onlyRetSpecs
          vars = {'Z1_Gaussians','Z2_Gaussians','Z3_Gaussians'};
          sizes = {[nIndep*M,1+2*kappa],[nIndep*M,1+2*kappa],[nIndep*M,1+kappa]};
          varargout = {vars,sizes};
          return;
       end

       % Generate numbers
       num = struct;
       num.Z1_Gaussians = randn(nIndep*M,1+2*kappa,precision);
       num.Z2_Gaussians = randn(nIndep*M,1+2*kappa,precision);
       num.Z3_Gaussians = randn(nIndep*M,1+kappa,precision);

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
           error(['MixedrBergomiClass:GenNumbersForMC: Pricer settings must ', ...
                  'be an object of type ''MonteCarloPricerSettingsClass''.']);
       end

        % Find the total number of steps needed:
        [~,nTotalSteps] = settings.GetNumSteps(ttm);
        
        % Generate random numbers:
        if ~isempty(seed);rng(seed);end
        settings.random = obj.GenNumbersForSim(settings.N,max(max(nTotalSteps)),...
                                               settings.price_estimation.antithetic,...
                                               settings.precision,false);

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
   %                                o V:          The instantaneous variance process.
   %                                o VIX2:       The VIX index squared.
   %                                o U1,U2       Let i = 1,2. We then approximate the kernel 
   %										      function K_i
   %
   %                                                K_i(t) approx(=) sum_{j=1}^m c_ij e^(-gamm_ij*t)
   %
   %                                              for t in [max(dt,kappa*dt),T] where 
   %
   %                                              T := max(T_basic,T_VIX2)
   %
   %                                              with 
   %
   %                                              T_basic = "maximum requested timepoint for 
   %                                                         variables S,V,U1,U2"
   %                                              T_VIX2 = "maximum requested timepoint for 
   %                                                        variable VIX2 plus 1/12" 
   %
   %                                              and where c = (c_i1,...,c_im), gamm = (gamm_i1,...
   %                                              gamm_im) are coefficients and dt the step size.
   %                                              The value of c, gamm (and m) will generally depend 
   %                                              on the particular i.
   %
   %                                              We then return the factors
   %
   %                                              U_ij(t) := 
   %                                                  eta_i*int_0^t exp(-gamm_ij*(t-s))dW_{i+1}(s)
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
   %     except for the parameter nI_VIX2 which is required if variable 'VIX2' is requested.
   %
   %   nI_VIX2:    [1x1 integer]             Number of integration points for computing the VIX
   %                                         index (squared). We here use a trapezoidal rule.
   %
   %   n_is:       [1x1 string]              Specifies how to interpret the 'n' input. Options are
   %                                         'per_unit_of_time' and 'total' (i.e. total number 
   %                                         of time points till the largest timepoint in the
   %                                         'timepoints' input). Default is 'per_unit_of_time'.
   %
   %   antithetic: [1x1 logical]             If true we also return antithetic sample paths where 
   %                                         the values of the Brownian vector W = (W1,W2,W3) are 
   %                                         flipped  to -W (still N paths in total). It is then  
   %                                         the first N/2 paths returned that will be the original 
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
   %                        If the simulation variable is 'Ui' for some i=1,2, we in addition 
   %                        include fields 'c' ([mx1 real])  'gamm' ([mx1 real]) for the 
   %                        coefficients explained under the 'variables' parameter.
   %
   % ---------------------------------------------------------------------------------------------

       % Parse inputs:
       inP = inputParser;
       addParameter(inP,'nI_VIX2',[]);
       addParameter(inP,'n_is','per_unit_of_time');
       addParameter(inP,'antithetic',false);
       addParameter(inP,'rndNumbers',[]);
       addParameter(inP,'tau',[]);
       parse(inP,varargin{:});
       v2struct(inP.Results);
       
       % Verify correlation matrix is ok:
       rho12 = obj.rho_12;rho13 = obj.rho_13;rho23 = obj.rho_23;
       if rho12^2 + rho13^2 + rho23^2 - 2*rho12*rho13*rho23 - eps > 1
           error('MixedrBergomiClass:Simulate: Invalid correlation matrix.');
       end
       
       % Check output variables are valid:
       expectedOutVars = {'S','V','VIX2','U1','U2'};
       if sum(ismember(variables,expectedOutVars)) ~= size(variables,2)
           error(['MixedrBergomiClass:Simulate: Some of the chosen output variables ',...
                  'are not supported.']);
       end
       
       % Check for duplicates:
       if size(unique(variables),2) ~= size(variables,2)
           error(['MixedrBergomiClass:Simulate: Values in ''variables'' parameter ',...
                  'must be unique.']);
       end
       
       if size(variables,1) ~= size(timepoints,1) || size(variables,2) ~= size(timepoints,2)
           error(['MixedrBergomiClass:Simulation: The inputs ''variables'' and ',...
                  '''timepoints'' must have the same dimensions']);
       end
       
       if antithetic && rem(N,2) ~= 0
           % Check divisibility by 2
           error(['MixedrBergomiClass:Simulate: When ',...
                  'antithetic samples are used we must have N divisible by 2.']);
       end       
       
       % Compute maximum time point:
       T = 0; % ~ horison to simulate till
       T_approx_K  = 0; % ~ horizon to approximate K's till
       for i=1:size(timepoints,2)
           if size(timepoints{i},1) > 1
               error(['MixedrBergomiClass:Simulate: The elements of ''timepoints'' ',...
                      'must all be row vectors.']);
           end
           if ~issorted(timepoints{i},'ascend')
               error(['MixedrBergomiClass:Simulate: Each element of ''timepoints'' ',...
                      'must be sorted in ascending order.']);
           end
           if any(timepoints{i} < 0)
               error(['MixedrBergomiClass:Simulate: The elements of ''timepoints'' ',...
                      'must be non-negative.']);
           end
           if size(timepoints{i},2) == 1 && timepoints{i} == 0
               error(['MixedrBergomiClass:Simulate: The elements of ''timepoints'' ',...
                      'are not allowed to be scalar [1x1 real] with a zero value.']);
           end
           if size(unique(timepoints{i}),2) ~= size(timepoints{i},2)
               error(['MixedrBergomiClass:Simulate: The elements of ''timepoints'' ',...
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
            error('MixedrBergomiClass:Simulate: Invalid value of input parameter ''n_is''.');
       end
       n_per_unit = 1/dt;
       grid_pts = dt*(0:nStepsTotalSim);
       T = max(grid_pts);
       
       % Unpack a few settings:
       simSettings = obj.pricerSettings.simulation;
       prec = obj.pricerSettings.precision;
       
       % Extract (rest of) model formulation:
       zeta = @(t)(obj.zeta.Eval(t));
       theta = @(t)(obj.theta.Eval(t));
       eta_1 = obj.eta_1;eta_2 = obj.eta_2;
       K1 = obj.K_1;K2 = obj.K_2;
       
       idxS = find(ismember(variables,'S'));
       idxV = find(ismember(variables,'V'));
       idxU1 = find(ismember(variables,'U1'));
       idxU2 = find(ismember(variables,'U2'));
       idxVIX2 = find(ismember(variables,'VIX2'));  
       
       if ~isempty(idxVIX2) && isempty(nI_VIX2)
           error(['MixedrBergomiClass:Simulate: When simulating ''VIX2'' we require ',...
                  'the input parameter ''nI_VIX2'' to be set.']);
       end
       
       % Compute sum of exponentials approximation:
       kernelApproxSettings = obj.pricerSettings.simulation.kernel_approximation;
       kappa = kernelApproxSettings.kappa;
       if kappa ~= 1 && ~isempty(idxVIX2)
           error(['MixedrBergomiClass:Simulate: We currently require ',...
                  'kappa = 1 when simulating VIX2.']);
       end    
       
       Npts = ceil2(n_per_unit*(T_approx_K-max(dt,kappa*dt))*kernelApproxSettings.n_multiplier);
       if Npts > kernelApproxSettings.n_cap
           Npts = kernelApproxSettings.n_cap;
           warning(['MixedrBergomiClass:Simulate: Number of sampling points for ',...
                    'sum of exponentials approximation was truncated at ', ...
                    num2str(kernelApproxSettings.n_cap), ' points.']);
       end       
       
       [c1,gamm1] = ApproximateKernel('BM2005','K',@(t)(K1.Eval(t)),'T',T_approx_K,...
                                    'delta',max(dt,dt*kappa),'N',Npts,...
                                    'epsilon',kernelApproxSettings.epsilon,...
                                    'error_measure',kernelApproxSettings.error_measure);
       [c2,gamm2] = ApproximateKernel('BM2005','K',@(t)(K2.Eval(t)),'T',T_approx_K,...
                                    'delta',max(dt,dt*kappa),'N',Npts,...
                                    'epsilon',kernelApproxSettings.epsilon,...
                                    'error_measure',kernelApproxSettings.error_measure);                               
       
       % Compute sum of exponentials approximation for K's squared:
       [c1_sqr,gamm1_sqr] = ApproximateKernel('BM2005','K',@(t)(K1.Eval(t).^2),'T',T_approx_K,...
                                    'delta',max(dt,dt*kappa),'N',Npts,...
                                    'epsilon',kernelApproxSettings.epsilon,...
                                    'error_measure',kernelApproxSettings.error_measure);
       [c2_sqr,gamm2_sqr] = ApproximateKernel('BM2005','K',@(t)(K2.Eval(t).^2),'T',T_approx_K,...
                                    'delta',max(dt,dt*kappa),'N',Npts,...
                                    'epsilon',kernelApproxSettings.epsilon,...
                                    'error_measure',kernelApproxSettings.error_measure);      
       
       % Determine time points to initially return U-values for:
      if ~any(ismember(variables,{'U1','U2','VIX2'}))
          % No time points needed:
          tU_sim = [];
          returnU = false;
      else
          % We return all time points needed for the requested variables:
          idxVolRelatedVars = find(ismember(variables,{'U1','U2','VIX2'}));
          tU_sim = [];
          for i=1:size(idxVolRelatedVars,2)
              tU_sim = unique([tU_sim,timepoints{idxVolRelatedVars(i)}]);
          end
          returnU = true;
      end 
       
       % Determine time points to initially return V-values for:
       if ~isempty(idxS)
           tV_sim = grid_pts;
       else
          idxVolVars = find(ismember(variables,{'V','U1','U2','VIX2'}));
          tV_sim = [];
          for i=1:size(idxVolVars,2)
              tV_sim = unique([tV_sim,timepoints{idxVolVars(i)}]);
          end
       end

       % Simulate random variables:
       if isempty(rndNumbers)
            rndNumbers = obj.GenNumbersForSim(N,nStepsTotalSim,antithetic,prec,false);
       end
       
       % Factorise correlation matrix:
       % Remark: W(t) = A*Z(t), AA^T = rho, A lower triangular.
       rho = [[1,rho12,rho13];...
              [rho12,1,rho23];...
              [rho13,rho23,1]];
       A = chol_mod(rho,10.^(-(16:-1:8)));
       
       % Make sure values are exactly zero in upper part of matrix so they can be detected as 
       % such by the 'SampleVolterra' function:
       A(1,2) = 0;A(1,3) = 0;A(2,3) = 0;
       
       % Sample Volterra integrals:
       Kid = KernelFunctionClass(1,0,@(obj,t)(1),0);
       W_bold = SampleVolterra(A,{Kid,K1,K2},[0,kappa,kappa],dt,antithetic,N,nStepsTotalSim,...
                                 {rndNumbers.numbers.Z1_Gaussians,...
                                  rndNumbers.numbers.Z2_Gaussians,...
                                  rndNumbers.numbers.Z3_Gaussians});
       
       % Simulate variance process:
       switch simSettings.scheme
           case 'hybrid_multifactor'
               % Remark: Input nStepsTotalSim/T may not be integer valued as the description of 
               % 'HybridMultifactorScheme' requires it to be. However, a careful reading of
               % the code shows that non-integer values are also valid.  
               [X1_nDrift,U1] = HybridMultifactorScheme(N,nStepsTotalSim/T,[],0,...
                                     0,eta_1,gamm1,c1,kappa,'K',K1,'returnU',returnU,...
                                     'positive',false,'tX',tV_sim,...
                                     'tU',tU_sim,'returndW',false,...
                                     'precision',prec,'explicit',simSettings.explicit,...
                                     'W_bold',W_bold{2});
               [X2_nDrift,U2] = HybridMultifactorScheme(N,nStepsTotalSim/T,[],0,...
                                          0,eta_2,gamm2,c2,kappa,...
                                          'K',K2,'returnU',returnU,...
                                          'positive',false,'tX',tV_sim,...
                                          'tU',tU_sim,'returndW',false,...
                                          'precision',prec,'explicit',simSettings.explicit,...
                                          'W_bold',W_bold{3});
               
           otherwise
               error('MixedrBergomiClass:Simulate: Scheme is not supported.');
       end
       
       % Compute and save variables:
       paths = struct;
       
       % V: 
       if ~isempty(idxV) || ~isempty(idxS) || ~isempty(idxVIX2)
           % Extract time points:
           tTruncV = X1_nDrift.t_truncated;
           tV = X1_nDrift.t;
           
           % Compute integrated K's (squared):           
           tSub = unique([0,tTruncV]');
           t_lb = tSub(1:end-1);
           t_ub = tSub(2:end);
           
           K1_squared_and_integrated_tmp = [0;cumsum(K1.Integrate(t_lb,t_ub,2,...
                                                          max(dt,dt*kappa),c1_sqr',gamm1_sqr'))]';
           K2_squared_and_integrated_tmp = [0;cumsum(K2.Integrate(t_lb,t_ub,2,...
                                                          max(dt,dt*kappa),c2_sqr',gamm2_sqr'))]';                                                      
               
           idxt = LocateTimepointsInDiscretisation([0;t_ub]',tTruncV);
           K1_squared_and_integrated = K1_squared_and_integrated_tmp(idxt);
           K2_squared_and_integrated = K2_squared_and_integrated_tmp(idxt);

           % Compute Y's:
           X_1 = X1_nDrift.values - 0.5*(eta_1^2)*K1_squared_and_integrated;
           X_2 = X2_nDrift.values - 0.5*(eta_2^2)*K2_squared_and_integrated;
           
           % Compute V:
           V = obj.c + zeta(tTruncV).*(theta(tTruncV).*exp(X_1) + (1-theta(tTruncV)).*exp(X_2));
           
           % Subset V:
           if ~isempty(idxV)
               idxt = LocateTimepointsInDiscretisation(tV,timepoints{idxV});
               paths.V.t = timepoints{idxV};
               paths.V.t_truncated = tTruncV(idxt);
               paths.V.values = V(:,idxt);
           end
       end
       
       % U's: 
       if ~isempty(idxU1)
           idxt = LocateTimepointsInDiscretisation(U1.t,timepoints{idxU1});
           paths.U1.t = timepoints{idxU1};
           paths.U1.t_truncated = U1.t_truncated(idxt);
           paths.U1.values = U1.values(:,:,idxt);
           paths.U1.c = c1;
           paths.U1.gamm = gamm1;
       end
       if ~isempty(idxU2)
           idxt = LocateTimepointsInDiscretisation(U2.t,timepoints{idxU2});
           paths.U2.t = timepoints{idxU2};
           paths.U2.t_truncated = U2.t_truncated(idxt);
           paths.U2.values = U2.values(:,:,idxt);
           paths.U2.c = c2;
           paths.U2.gamm = gamm2;
       end    
              
       % S: 
       if ~isempty(idxS)
           dlogS = zeros(N,size(V,2),prec);
           dlogS(:,2:end) = sqrt(V(:, 1:end-1)).*W_bold{1}(:,1:size(V,2)-1) ...
                            - 0.5*V(:,1:end-1)*dt; 
           S = exp(cumsum(dlogS,2));
           idxt = LocateTimepointsInDiscretisation(tV,timepoints{idxS});
           paths.S.t = timepoints{idxS};
           paths.S.t_truncated = tTruncV(idxt);
           paths.S.values = S(:,idxt);
       end
       
       % Remark: For the code below we currently require kappa = 1. That indeed kappa = 1, is 
       % verified further up in the code.       

       % VIX2:
       if ~isempty(idxVIX2)
           % Compute time points:
            idxtU = LocateTimepointsInDiscretisation(U1.t,timepoints{idxVIX2});
            idxtV = LocateTimepointsInDiscretisation(tV,timepoints{idxVIX2});
            paths.VIX2.t = timepoints{idxVIX2};
            paths.VIX2.t_truncated = U1.t_truncated(idxtU);        
            paths.VIX2.values = zeros(N,size(paths.VIX2.t,2));
            
            % Compute integration points (i.e. time to expiry of xi's):
            DELTA = 1/12;
            dTau = DELTA/(nI_VIX2-1);
            tauIntPts = (0:dTau:DELTA);
            idxLow = tauIntPts <= dt;
            idxHigh = ~idxLow;
            tauPtsLow = tauIntPts(idxLow);
            tauPtsHigh = tauIntPts(idxHigh);
            
            % We need to know xi at maturity 'dt' for linear interpolation in front end.
            % Thus include in vector:
            tauPtsHighAdd = [dt,tauPtsHigh]; 
            nPtsLow = size(tauPtsLow,2);
           
            % Compute variables shared by all time points:
            dummy1 = c1.*exp(-gamm1*tauPtsHighAdd);
            dummy2 = c2.*exp(-gamm2*tauPtsHighAdd);
            
            xiVals = zeros(N,size(tauIntPts,2));
            for i=1:size(paths.VIX2.t,2)
                f1 = 0.5*(eta_1^2)*K1.Integrate(tauPtsHighAdd',...
                                            paths.VIX2.t(i) + tauPtsHighAdd',...
                                            2,max(dt,dt*kappa),c1_sqr',gamm1_sqr').';
                f2 = 0.5*(eta_2^2)*K2.Integrate(tauPtsHighAdd',...
                                            paths.VIX2.t(i) + tauPtsHighAdd',...
                                            2,max(dt,dt*kappa),c2_sqr',gamm2_sqr').';                                          

                intVals_1 = U1.values(:,:,idxtU(i))*dummy1 - f1;
                intVals_2 = U2.values(:,:,idxtU(i))*dummy2 - f2;
                
                % Compute forward variances at longer maturities:
                zetaVals = zeta(paths.VIX2.t(i) + tauPtsHigh);
                thetaVals = theta(paths.VIX2.t(i) + tauPtsHigh);
                xiVals(:,nPtsLow+1:end) = obj.c + zetaVals.*(thetaVals.*exp(intVals_1(:,2:end)) ...
                                                       + (1-thetaVals).*exp(intVals_2(:,2:end)));
                
                % Compute xi(dt) value:
                t_plus_dt = paths.VIX2.t(i)+dt;
                xi_dt = obj.c + zeta(t_plus_dt).*(theta(t_plus_dt).*exp(intVals_1(:,1)) ...
                              + (1-theta(t_plus_dt)).*exp(intVals_2(:,1)));
                
                % Compute forward variances at shorter maturities:
                xiVals(:,1:nPtsLow) = V(:,idxtV(i)).*((dt-tauPtsLow)./dt) + xi_dt.*(tauPtsLow./dt);
                
                % Integrate over forward variances:                                   
                paths.VIX2.values(:,i) = (1/DELTA)*(0.5*xiVals(:,1) + sum(xiVals(:,2:end-1),2)...
                                                    + 0.5*xiVals(:,end))*dTau;
                                                
            end

       end

   end
end

end

