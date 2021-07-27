classdef MixedQuadraticRoughHyperbolicModelClass < PricingModelClass & handle
%
%   Implements the mixed quadratic rough hyperbolic model. Let r(t) and delta(t) respectively denote 
%   the risk-free interest rate and continuous dividends yield. Both are assumed deterministic 
%   functions of time. The asset price S(t) is assumed to have risk-neutral dynamics of the form
%
%     dS(t) = S(t)(r(t) - delta(t))dt + S(t)sqrt(V(t))dW_1(t)
%
%   where
%
%     V(t) = nu(t)*(theta*X1(t) + (1-theta)*X2(t)) + c,
%
%   with
%
%     f(y) = y + sqrt(y^2 + 1),
%
%     X1(t) = f(Y1(t)^2 - b1) - f(-b1)
%     X2(t) = f(Y2(t)^2 - b2) - f(-b2)
%
%     Y1(t) = zeta1(t) + eta1*int_0^t K1(t-s) dW_2(s)
%     Y2(t) = zeta2(t) + eta2*int_0^t K2(t-s) dW_3(s)
%
%   Here zeta1, zeta2, K1 and K2 are deterministic functions (K's completely monotone). 
%   Parameters eta1, eta2, c are non-negative and theta must lie in [0,1] (no restrictions apply 
%   to b1 and b2). The Brownian vector W = (W1,W2,W3) has a general correlation structure as 
%   specified by the parameters rho_12, rho_13, rho_23. 
% 
% -----------------------------------------------------------------------------------------------
%  Properties  
% -----------------------------------------------------------------------------------------------
%   theta:    [1x1 real]                    Mixing weight.
%   c:        [1x1 real]                    Minimum instantaneous variance.
%   b1:       [1x1 real]                    Shift parameter for X1.
%   b2:       [1x1 real]                    Shift parameter for X2.
%   eta1:     [1x1 real]                    Volatility of Y1.
%   eta2:     [1x1 real]                    Volatility of Y2.
%   rho_12:   [1x1 real]                    Correlation between W1 and W2.
%   rho_13:   [1x1 real]                    Correlation between W1 and W3.
%   rho_23:   [1x1 real]                    Correlation between W2 and W3.
%   K1:       [1x1 KernelFunctionClass]     Volterra kernel for Y1.
%   K2:       [1x1 KernelFunctionClass]     Volterra kernel for Y2.
%   zeta1:    [1x1 CurveClass]              Term structure of Y1. 
%   zeta2:    [1x1 CurveClass]              Term structure of Y2. 
%	nu:		  [1x1 CurveClass]				Term structure of overall process.
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
%         'asset_price'.
%       o obj.pricerSettings.price_estimation.VIX.conditional_monte_carlo: Must be set to false.
%       o obj.pricerSettings.price_estimation.VIX.control_variate: Supported values are 'none'.
%
% -----------------------------------------------------------------------------------------------

properties
    theta
    c
    b1
    b2
    eta1
    eta2
    rho_12
    rho_13
    rho_23
    K1
    K2
    zeta1
    zeta2
	nu
end

methods
   function obj = MixedQuadraticRoughHyperbolicModelClass(varargin)
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
            
            sim.Ehyp.dx = 0.01;
            sim.Ehyp.interpolation = 'linear';
            
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

        % Set or adjust the input curves:
        if isprop(obj,'zeta1') && isempty(obj.zeta1)
            obj.zeta1 = CurveClass('gridpoints',1,'values',0);
        end
        if isprop(obj,'zeta1') && isnumeric(obj.zeta1)
            obj.zeta1 = CurveClass('gridpoints',1,'values',obj.zeta1);
        end
        
        if isprop(obj,'zeta2') && isempty(obj.zeta2)
            obj.zeta2 = CurveClass('gridpoints',1,'values',0);
        end
        if isprop(obj,'zeta2') && isnumeric(obj.zeta2)
            obj.zeta2 = CurveClass('gridpoints',1,'values',obj.zeta2);
        end  
		
        if isprop(obj,'nu') && isempty(obj.nu)
            obj.nu = CurveClass('gridpoints',1,'values',1);
        end
        if isprop(obj,'nu') && isnumeric(obj.nu)
            obj.nu = CurveClass('gridpoints',1,'values',obj.nu);
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

       
       if condMC_S
           error(['MixedQuadraticRoughHyperbolicModelClass:GetPricesSub: ',...
                  'Conditional Monte Carlo is not supported for options on S(t).']);
       end       
       if condMC_V
           error(['MixedQuadraticRoughHyperbolicModelClass:GetPricesSub: ',...
                  'Conditional Monte Carlo is not supported for VIX options.']);
       end
       if ~(strcmpi(cvS,'none') || strcmpi(cvS,'asset_price'))
           error(['MixedQuadraticRoughHyperbolicModelClass:GetPricesSub: ',...
                  'obj.price_estimation.S.control_variate only supports options ''none''',...
                  'and ''asset_price''.']);
       end
       if ~strcmpi(cvV,'none')
           error(['MixedQuadraticRoughHyperbolicModelClass:GetPricesSub: ',...
                  'obj.price_estimation.VIX.control_variate ',...
                  'only supports option ''none']);
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
       [sim_variables, tOut] = deal({});
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
                cvMCEst = cvV;

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
           error(['MixedQuadraticRoughHyperbolicModelClass:GetNumbersForSim: Scheme ', ...
               scheme, ' is not supported.']);
       end
       
       kappa = obj.pricerSettings.simulation.kernel_approximation.kappa;

       if anti
           % Check divisibility by 2
           if rem(N,2) ~= 0
               error(['MixedQuadraticRoughHyperbolicModelClass:GenNumbersForSim: When ',...
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
           error(['MixedQuadraticRoughHyperbolicModelClass:GenNumbersForMC: Pricer settings ', ...
                  'must be an object of type ''MonteCarloPricerSettingsClass''.']);
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
   %                                o X1,X2:      See main description.
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
   %                                         the values of the W(t) Brownian motion are flipped 
   %                                         to -W(t) (still N paths in total). It is then the 
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
           error('MixedQuadraticRoughHyperbolicModelClass:Simulate: Invalid correlation matrix.');
       end
       
       % Check output variables are valid:
       expectedOutVars = {'S','V','VIX2','X1','X2'};
       if sum(ismember(variables,expectedOutVars)) ~= size(variables,2)
           error(['MixedQuadraticRoughHyperbolicModelClass:Simulate: ',...
                  'Some of the chosen output variables ',...
                  'are not supported.']);
       end
       
       % Check for duplicates:
       if size(unique(variables),2) ~= size(variables,2)
           error(['MixedQuadraticRoughHyperbolicModelClass:Simulate: Values in ',...
                  '''variables'' parameter must be unique.']);
       end
       
       if size(variables,1) ~= size(timepoints,1) || size(variables,2) ~= size(timepoints,2)
          error(['MixedQuadraticRoughHyperbolicModelClass:Simulation: The inputs ''variables'' ',...
                 'and ''timepoints'' must have the same dimensions']);
       end
       
       if antithetic
           % Check divisibility by 2
           if rem(N,2) ~= 0
               error(['MixedQuadraticRoughHyperbolicModelClass:Simulate: When ',...
                      'antithetic samples are used we must have N divisible by 2.']);
           end
           Nindep = N/2;
       else
           Nindep = N;
       end              
       
       % Compute maximum time point:
       T = 0; % ~ horison to simulate till
       T_approx_K  = 0; % ~ horizon to approximate K's till
       for i=1:size(timepoints,2)
           if size(timepoints{i},1) > 1
               error(['MixedQuadraticRoughHyperbolicModelClass:Simulate: The elements of ',...
                      '''timepoints'' must all be row vectors.']);
           end
           if ~issorted(timepoints{i},'ascend')
               error(['MixedQuadraticRoughHyperbolicModelClass:Simulate: Each element of ',...
                      '''timepoints'' must be sorted in ascending order.']);
           end
           if any(timepoints{i} < 0)
               error(['MixedQuadraticRoughHyperbolicModelClass:Simulate: The elements of ',...
                      '''timepoints'' must be non-negative.']);
           end
           if size(timepoints{i},2) == 1 && timepoints{i} == 0
               error(['MixedQuadraticRoughHyperbolicModelClass:Simulate: The elements of ',...
                      '''timepoints'' are not allowed to be scalar [1x1 real] with a zero value.']);
           end
           if size(unique(timepoints{i}),2) ~= size(timepoints{i},2)
               error(['MixedQuadraticRoughHyperbolicModelClass:Simulate: The elements of ',...
                      '''timepoints'' must contain distinct values only.']);               
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
            error(['MixedQuadraticRoughHyperbolicModelClass:Simulate: Invalid value of input',...
                  'parameter ''n_is''.']);
       end
       n_per_unit = 1/dt;
       grid_pts = dt*(0:nStepsTotalSim);
       T = max(grid_pts);
       
       % Unpack a few settings:
       simSettings = obj.pricerSettings.simulation;
       prec = obj.pricerSettings.precision;
       
       % Extract model formulation:
       zeta1 = @(t)(obj.zeta1.Eval(t));
       zeta2 = @(t)(obj.zeta2.Eval(t));
	   nu = @(t)(obj.nu.Eval(t));
       theta = obj.theta;
       eta1 = obj.eta1;eta2 = obj.eta2;       
       b1 = obj.b1;b2 = obj.b2;
       K1 = obj.K1;K2 = obj.K2;
       
       idxS = find(ismember(variables,'S'));
       idxV = find(ismember(variables,'V'));
       idxVIX2 = find(ismember(variables,'VIX2'));
       idxX1 = find(ismember(variables,'X1'));
       idxX2 = find(ismember(variables,'X2'));
       
       if ~isempty(idxVIX2) && isempty(nI_VIX2)
           error(['MixedQuadraticRoughHyperbolicModelClass:Simulate: When simulating ''VIX2'' ',...
                  'we require the input parameter ''nI_VIX2'' to be set.']);
       end
       
       % Compute sum of exponentials approximation:
       kernelApproxSettings = obj.pricerSettings.simulation.kernel_approximation;
       kappa = kernelApproxSettings.kappa;
       if kappa ~= 1 && ~isempty(idxVIX2)
           error(['MixedQuadraticRoughHyperbolicModelClass:Simulate: We currently require ',...
                  'kappa = 1 when simulating VIX2.']);
       end    
       
       Npts = ceil2(n_per_unit*(T_approx_K-max(dt,kappa*dt))*kernelApproxSettings.n_multiplier);
       if Npts > kernelApproxSettings.n_cap
           Npts = kernelApproxSettings.n_cap;
           warning(['MixedQuadraticRoughHyperbolicModelClass:Simulate: Number of sampling ',...
                    'points for sum of exponentials approximation was truncated at ', ...
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

       % Determine time points to initially return U-values for:
      if ~any(ismember(variables,{'VIX2'}))
          % No time points needed:
          tU_sim = [];
          returnU = false;
      else
          % We return all time points needed for the requested variables:
          idxVolRelatedVars = find(ismember(variables,{'VIX2'}));
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
          idxVolVars = find(ismember(variables,{'V','VIX2'}));
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
               [Y1paths,U1paths] = HybridMultifactorScheme(N,nStepsTotalSim/T,[],zeta1,...
                                        0,eta1,gamm1,c1,kappa,'K',K1,'returnU',returnU,...
                                        'positive',false,'tX',tV_sim,...
                                        'tU',tU_sim,'returndW',false,...
                                        'precision',prec,'explicit',simSettings.explicit,...
                                        'W_bold',W_bold{2});

               U1paths.c = c1;
               U1paths.gamm = gamm1;                                    
                                    
               [Y2paths,U2paths] = HybridMultifactorScheme(N,nStepsTotalSim/T,[],zeta2,...
                                        0,eta2,gamm2,c2,kappa,'K',K2,'returnU',returnU,...
                                        'positive',false,'tX',tV_sim,...
                                        'tU',tU_sim,'returndW',false,...
                                        'precision',prec,'explicit',simSettings.explicit,...
                                        'W_bold',W_bold{3});

               U2paths.c = c2;
               U2paths.gamm = gamm2;
                              
           otherwise
               error('MixedQuadraticRoughHyperbolicModelClass:Simulate: Scheme is not supported.');
       end
       
       % Compute and save variables:
       paths = struct;
       
       % V: 
       if ~isempty(idxV) || ~isempty(idxS) || ~isempty(idxVIX2)
           % Compute V:
           Vpaths = struct;
           Vpaths.t = Y1paths.t;
           Vpaths.t_truncated = Y1paths.t_truncated;
           fhyp = @(y)(y + sqrt(y.^2 + 1));
           X1 = fhyp( Y1paths.values.^2 - b1) - fhyp(-b1);
           X2 = fhyp( Y2paths.values.^2 - b2) - fhyp(-b2);
           Vpaths.values = max(nu(Y1paths.t_truncated).*(theta*X1 + (1-theta)*X2) + obj.c,0); 
		   % Remark for the above: max() is only there to handle problems with round-off errors.

           % Extract selected timepoints:
           if ~isempty(idxV)
               idxt = LocateTimepointsInDiscretisation(Vpaths.t,timepoints{idxV});
               paths.V.t = timepoints{idxV};
               paths.V.t_truncated = Vpaths.t_truncated(idxt);
               paths.V.values = Vpaths.values(:,idxt);
           end
           
       end
       
       if ~isempty(idxX1)
           idxt = LocateTimepointsInDiscretisation(Vpaths.t,timepoints{idxX1});
           paths.X1.t = timepoints{idxX1};
           paths.X1.t_truncated = Vpaths.t_truncated(idxt);
           paths.X1.values = X1(:,idxt);
       end
       
       if ~isempty(idxX2)
           idxt = LocateTimepointsInDiscretisation(Vpaths.t,timepoints{idxX2});
           paths.X2.t = timepoints{idxX2};
           paths.X2.t_truncated = Vpaths.t_truncated(idxt);
           paths.X2.values = X2(:,idxt);           
       end       

       % S: 
       if ~isempty(idxS)             
           dlogS = zeros(N,size(Vpaths.values,2),prec);
           dlogS(:,2:end) = sqrt(Vpaths.values(:, 1:end-1))...
                                  .*W_bold{1}(:,1:size(Vpaths.values,2)-1) ...
                                  - 0.5 *Vpaths.values(:,1:end-1)*dt; 
           S = exp(cumsum(dlogS,2));
           idxt = LocateTimepointsInDiscretisation(Vpaths.t,timepoints{idxS});
           paths.S.t = timepoints{idxS};
           paths.S.t_truncated = Vpaths.t_truncated(idxt);
           paths.S.values = S(:,idxt);
       end
       
       % Remark: For the code below we currently require kappa = 1. That indeed kappa = 1, is 
       % verified further up in the code. 

       % VIX2:
       if ~isempty(idxVIX2)
           % Compute time points:
            idxtU = LocateTimepointsInDiscretisation(U1paths.t,timepoints{idxVIX2});
            idxtV = LocateTimepointsInDiscretisation(Vpaths.t,timepoints{idxVIX2});
            paths.VIX2.t = timepoints{idxVIX2};
            paths.VIX2.t_truncated = U1paths.t_truncated(idxtU);        
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
            
            % Compute integrated K's (squared):
            [K1_sqr_int, K2_sqr_int] = deal(zeros(size(tauPtsHighAdd),prec));
            for i=1:numel(tauPtsHighAdd)
                K1_sqr_int(i) = obj.K1.IntegrateVolterra(0,tauPtsHighAdd(i),...
                                                         tauPtsHighAdd(i),2,@(s)(1));
                K2_sqr_int(i) = obj.K2.IntegrateVolterra(0,tauPtsHighAdd(i),...
                                                         tauPtsHighAdd(i),2,@(s)(1));
            end
           
            % Compute variables shared by all time points:
            dummy1 = c1.*exp(-gamm1*tauPtsHighAdd);
            dummy2 = c2.*exp(-gamm2*tauPtsHighAdd);
            
            b_tilde_1 = eta1*sqrt(K1_sqr_int);
            b_tilde_2 = eta2*sqrt(K2_sqr_int);
            c_tilde_1 = -b1;
            c_tilde_2 = -b2;

            % Settings for hyperbolic expectation:
            dx = obj.pricerSettings.simulation.Ehyp.dx;
            imethod = obj.pricerSettings.simulation.Ehyp.interpolation;
            
            % Compute all a-tilde values (+ ranges at all time-to-maturities):
            [a_tilde_1,a_tilde_2] = deal(cell(size(paths.VIX2.t,2),numel(tauPtsHighAdd)));
            [a_tilde_min,a_tilde_max] = deal(zeros(numel(tauPtsHighAdd),2));
            a_tilde_min(:) = Inf;a_tilde_max(:) = -Inf;
            for i=1:size(paths.VIX2.t,2)
                % Compute forward variances at longer maturities:
                eta_times_stochInt_1 = U1paths.values(:,:,idxtU(i))*dummy1;
                eta_times_stochInt_2 = U2paths.values(:,:,idxtU(i))*dummy2;
                for j=1:numel(tauPtsHighAdd)
                    % Compute a tilde:
                    a_tilde_1{i,j} = zeta1(paths.VIX2.t(i) + tauPtsHighAdd(j)) ...
                                                           + eta_times_stochInt_1(:,j);
                    a_tilde_2{i,j} = zeta2(paths.VIX2.t(i) + tauPtsHighAdd(j)) ...
                                                           + eta_times_stochInt_2(:,j);
                    a_tilde_min(j,1) = min([a_tilde_min(j,1),min(min(a_tilde_1{i,j}))]);
                    a_tilde_min(j,2) = min([a_tilde_min(j,2),min(min(a_tilde_2{i,j}))]);
                    a_tilde_max(j,1) = max([a_tilde_max(j,1),max(max(a_tilde_1{i,j}))]);
                    a_tilde_max(j,2) = max([a_tilde_max(j,2),max(max(a_tilde_2{i,j}))]);
                end
            end
            
            % Compute interpolation functions for hyperbolic expectation at all time-to-maturities:
            [a_tilde_val,Ehyp_val] = deal(cell(numel(tauPtsHighAdd),2));
            for i=1:numel(tauPtsHighAdd)
              a_tilde_val{i,1} = unique([(a_tilde_min(i,1):dx:a_tilde_max(i,1))';a_tilde_max(i,1)]);
              intFun = @(x)((fhyp((a_tilde_val{i,1} + b_tilde_1(i).*x).^2 + c_tilde_1)...
                                   ./sqrt(2*pi)).*exp(-x.^2/2));
              Ehyp_val{i,1} = integral(intFun,-Inf,Inf,'ArrayValued',true);
              a_tilde_val{i,2} = unique([(a_tilde_min(i,2):dx:a_tilde_max(i,2))';a_tilde_max(i,2)]);
              intFun = @(x)((fhyp((a_tilde_val{i,2} + b_tilde_2(i).*x).^2 + c_tilde_2)...
                                   ./sqrt(2*pi)).*exp(-x.^2/2));
              Ehyp_val{i,2} = integral(intFun,-Inf,Inf,'ArrayValued',true);                
            end
            
            % Compute actual forward variances and VIX values:
            xiVals = zeros(N,size(tauIntPts,2));
            for i=1:size(paths.VIX2.t,2)
                % Compute forward variances at longer maturities:
                for j=1:numel(tauPtsHighAdd)
                    Ehyp_val_1 = interp1(a_tilde_val{j,1},Ehyp_val{j,1},a_tilde_1{i,j},imethod);
                    Ehyp_val_2 = interp1(a_tilde_val{j,2},Ehyp_val{j,2},a_tilde_2{i,j},imethod);
                    
					% Compute nu() at relevant time point:
					nu_val = nu(paths.VIX2.t(i) + tauPtsHighAdd(j));
					
                    % Compute expectation:
                    if j == 1
                        xi_dt = nu_val.*(...
                                  obj.theta*(Ehyp_val_1 - fhyp(-obj.b1))...
                             + (1-obj.theta)*(Ehyp_val_2 - fhyp(-obj.b2)))...
                             + obj.c;
                    else
                        xiVals(:,nPtsLow+j-1) = nu_val.*(...
                                 obj.theta*(Ehyp_val_1 - fhyp(-obj.b1))...
                              + (1-obj.theta)*(Ehyp_val_2 - fhyp(-obj.b2))...
                              ) + obj.c;
                    end
					                    
                end

                % Compute forward variances at shorter maturities:
                xiVals(:,1:nPtsLow) = Vpaths.values(:,idxtV(i)).*((dt-tauPtsLow)./dt)...
                                                                   + xi_dt.*(tauPtsLow./dt);                
                                  
                                                              
                % Integrate over forward variances:                                   
                paths.VIX2.values(:,i) = max((1/DELTA)*(0.5*xiVals(:,1) ...
                                                    + sum(xiVals(:,2:end-1),2)...
                                                    + 0.5*xiVals(:,end))*dTau,0);  
                                                
            end

       end

   end

        
end

end

