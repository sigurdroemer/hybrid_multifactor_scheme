classdef MonteCarloPricerSettingsClass < GenericSuperClass
%
% 	Class to hold settings for a Monte Carlo based pricer as well as some useful functionality.
%
% -----------------------------------------------------------------------------------------------
%  Properties
% -----------------------------------------------------------------------------------------------
%  tn,n:    [Mx1 real resp.   Properties obj.tn and obj.n joint specify the desired number of steps 
%             ...  integer]   per year for different intervals of expiries. The first interval of 
%                             expiries goes from 0 to obj.tn(1). The second from obj.tn(1) to 
%                             obj.tn(2), etc. The desired number of steps per year for each group 
%                             is then specified by the property obj.n. 
%
%                             Remarks:
%                               o When Richardson-Romberg (RR) extrapolation is NOT used, obj.n 
%                                 specifies the actual number of steps per year. Let T be the 
%                                 vector of expiries for a particular expiry group (i.e. all 
%                                 corresponding to a particular interval as specified in obj.tn). 
%                                 Let n be the corresponding number of steps per year as given in 
%                                 obj.n. The step size is then dt = 1/n and the simulation grid
%                                 thus given by
%
%                                   0 < 1/n < 2/n < ... < floor(max(T)*n)/n
%
%                                 Expiries are then priced by using simulated values at the nearest 
%                                 grid point just below.
%
%                               o When RR extrapolation IS used we need simulated values under both 
%                                 a fine and coarse discretisation. In this case obj.n specifies 
%                                 the desired number of steps per year for the finer of the two 
%                                 discretisations. Also, the actual number of steps per year (for 
%                                 the finer discretisation) may be slightly higher than what is 
%                                 specified in obj.n. The reason is that we need the total number 
%                                 of steps till the last expiry in a given expiry-group to be 
%                                 divisible by 2. Letting dt <= 1/n be the actual step size used, 
%                                 we use the simulation grid
%
%                                    0 < dt < 2*dt < ...
%
%                                 and once again price expiries by using simulated values at the 
%                                 nearest grid point just below.
%
%                                 It is briefly explained how 'dt' is computed: Defining
%
%                                    total_num_steps_initial := floor(max(T)*n)      
%
%                                 we adjust that number to
%
%                                    total_num_steps = 2*ceil(total_num_steps_initial/2)
%
%                                 giving us the actual (total) number of steps to be used for the 
%                                 finer discretisation. The step size is then computed as
%
%                                    dt = max(T)/total_num_steps                            (1)
% 
%                                 for the finer discretisation, and as
%
%                                    dt = max(T)/(total_num_steps/2)                        (2)
%
%                                 for the coarser discretisation.                   
%
%                                 See the explanation under the property 'richardson_romberg' for 
%                                 more information on how RR extrapolation is used (see field 
%                                 'price_estimation' further down).
%   
%                               o Code for computing the actual number of steps is contained 
%                                 in method 'GetNumSteps'.
%
%  N:          [1x1 integer]  Assuming numRepeat = 1 (see below) property 'N' sets the number of 
%                             paths to be used (including any antithetic ones). When numRepeat > 1 
%                             the total number of paths is N*numRepeat as we run the pricing 
%                             algorithm a 'numRepeat' number of times each time using 'N' paths.
%
%  numRepeat:  [1x1 integer]  Number of times to repeat the computations. Prices obtained from each 
%                             run are then averaged. Total number of paths used is therefore 
%                             N*numRepeat. Useful to limit memory usage. When control variates are 
%                             used, it is though highly recommended to not keep N too low to avoid 
%                             biasing results.
% 
%  o seed:     [1x1 integer   If non-empty we call rng(seed) before generating random numbers 
%               ... or empty] with 'GenNumbersForMC'. Not supported if numRepeat > 1.
%
%  simulation: [1x1 struct]   Struct containing settings on how paths should be simulated. The 
%                             following fields apply:
%
%                             scheme:  [1x1 string] Options are 'locally_lognormal', 
%                                                   'hybrid_multifactor'.
%
%                             epsilon: [1x1 real]   Applies to scheme 'locally_lognormal'. See 
%                                                   function 'LocallyLognormalScheme'. Can be left 
%                                                   empty in which case the default value specified
%                                                   in that function is used.
%
%                             nu:      [1x1 real]   Applies to scheme 'locally_lognormal'. See 
%                                                   function 'LocallyLognormalScheme'. Can be left 
%                                                   empty in which case the default value specified
%                                                   in that function is used.
%
%                             tau_max: [1x1 real]   Applies to scheme 'locally_lognormal'. See 
%                                                   function 'LocallyLognormalScheme'. Can be left 
%                                                   empty in which case the default value specified
%                                                   in that function is used.
%
%                             kernel_approximation: [1x1 struct] 
%
%                                 Contains settings for numerically approximating the kernel 
%                                 function(s).
%
%                                 Letting dt > 0 be the simulation step size and kappa a 
%                                 non-negative integer, each kernel is approximated exactly on the 
%                                 interval [0,kappa*dt] and by a sum of exponential functions
%                                 on the domain from max(dt,kappa*dt) till the relevant upper end  
%                                 point. The 'BM2005' method of the function 'ApproximateKernel' is 
%                                 used to determine the sum of exponentials approximation. 
%
%                                 Remark: The lower end point for the sum of exponentials
%                                 approximation is computed as max(dt,dt*kappa) since the 'BM2005'
%                                 method does not work when the domain includes a singularity 
%                                 (possibly the case when kappa = 0).
%
%                                 The following subfields apply and determine how the approximation
%                                 is computed.
%       
%                                    kappa:         [1x1 integer]
%                                    n_multiplier:  [1x1 integer] 
%                                    n_cap:         [1x1 integer]
%                                    error_measure: [1x1 string]
%                                    epsilon:       [1x1 real]
%
%                                 The 'kappa' field is self-explanatory.
%
%                                 Parameters 'n_multiplier' and 'n_cap' sets the number of points
%                                 to be used when running the method 'BM2005'. More precisely, 
%                                 say that the simulation algorithm wishes to approximate K
%                                 on an interval [a,b] and that the simulation grid uses n_per_unit
%                                 number of steps per unit of time (possibly non-integer). We then 
%                                 set the total number of sampling points to be used for the BM2005 
%                                 method as
%
%                                   num. points = min(ceil2(n_per_unit*(b-a)*n_multiplier), n_cap)
%
%                                 For fields 'error_measure' and 'epsilon', see function 
%                                 'ApproximateKernel' for an explanation.
%
%                            Additional restrictions may apply for a given model; see the 
%                            particular class description.
%
%  price_estimation: [1x1 struct] Settings for price estimation. Must contain the following fields:
%
%      o antithetic:         [1x1 logical]  
%           
%             If true half of simulated paths will be antithetic.
%
%      o richardson_romberg: [1x1 logical]  
%           
%             If true we use Richardson-Romberg (RR) extrapolation. Let X(t) be a placeholder 
%             stochastic process representing either S(t) or VIX(t) and let f be the payoff 
%             function of either a put or call option. 
%
%             Say we wish to price an option with expiry Ti > 0 and payoff function f. We then 
%             hypothesize the relationship                                                                     
%
%                  E(f(X(Ti)^m)) = E(f(X(Ti)) + cm^(-p1)  + o(m^(-p2))                      (3)
%
%             where c is some constant, 0 < p1 < p2, and X(Ti)^m denotes the numerical solution of 
%             X(Ti) using m steps in total till the expiry Ti.
%
%             To apply RR extrapolation, we sample X using both some m1 number of steps and some 
%             m2 number of steps. Here m1 < m2. Ideally m2 = 2*m1 but that is only approximately
%             true in practise.
%
%             Recalling equations (1) and (2) from the remarks under the description of 
%             properties 'tn' and 'n', we more precisely have
%
%                  m2 = floor(Ti/dt)
%
%             where dt refers to the value from equation (1), and
%   
%                  m1 = floor(Ti/dt)
%
%             where dt instead refers to the value from (2). Note that m1 and m2 both depend on 
%             'n' (the desired number of steps per year for the particular group of expiries). 
%             We may therefore write m1(n) and m2(n). See also that 
%
%                  m1(n),m2(n) -> infinity                                                  (4)
%
%             and
%
%                  m2(n)/m1(n) -> 2                                                         (5)
%
%             both as n -> infinity.
%   
%             Using X(Ti)^m1 and X(Ti)^m2 we construct an extrapolated estimator as 
%
%                  (1/(a-1))*( a*E(f(X(Ti)^m2)) - E(f(X(Ti)^m1))) = E(f(X(Ti))) + o(m2(n)^(-p2))
%
%             where we have defined a := (m2/m1)^p1. The equality in the above follow due to our 
%             hypothesis (3) as well as the limits (4)-(5). Note that the error o(m2(n)^(-p2)) may 
%             be replaced by o(m1(n)^(-p2)) due to the limit (5).
%
%             Under the hypothesis, the rate of convergence of the extrapolated estimator is 
%             therefore p2 (> p1).
%
%             In practise, the expectations operator E(...) is of course replaced by its empirical
%             version.
%
%             We do not use RR extrapolation to obtain the VIX futures prices. In this case we 
%             simply use the estimates obtained by using the finer of the two discretisations.
%
%             If RR extrapolation is used we require the field 'weak_error_rate' under 'S' and 
%             'VIX' to be specified.
%
%             	o S:          [1x1 struct]   Additional settings for pricing options on the 
%                                            underlying asset.
%
%               o VIX:        [1x1 struct]   Additional settings for pricing VIX options.
%
%              For the last two fields 'S' and 'VIX' the following subfields apply:
%
%                 o control_variate:         [1x1 string]   Specifies any control variates to use.  
%                                                           For field 'S' options are: 'asset_price' 
%                                                           ,'timer_option', 'none'. For field 'VIX'
%                                                           options are 'VIX2', 'none'.
%
%                 o conditional_monte_carlo: [1x1 logical]  Set to true to use conditional Monte 
%                                                           Carlo.
%
%                 o weak_error_rate:         [1x1 function] Function returning value of 'p1' as 
%                                                           explained under 'richardson_romberg'. 
%                                                           Function must take a PricingModelClass 
%                                                           object as input and then return a 
%                                                           [1x1 real] value giving 'p1'.
%
%                                                           This property is only required if RR 
%                                                           extrapolation is used.
%
%                 o option_type:             [1x1 string]   Specifies which option type should be 
%                                                           used to estimate prices. Possible 
%                                                           choices are 'call', 'put' and 'otm' 
%                                                           (out-of-the-money) options. The latter 
%                                                           is recommended.
%
%     Additional restrictions may apply to the 'price_estimation' field; consult the class 
%     description of the particular model.
%
%  o precision:   [1x1 string]   Computational precision. Options are 'single' and 'double'.
%
%  o random:      [(depends)]    Random numbers for simulations. Allows reuse if so desired. 
%                                Not intended for direct use by an end-user but may instead
%                                be used internally in code.
%
%  o optim:       [(depends)]    Stores temporary information when calibrating model.
%
% -----------------------------------------------------------------------------------------------
    
properties
    tn
    n
    n_VIX2_integrate
    n_EVIX2_integrate
    N
    numRepeat
    seed
    simulation
    price_estimation
    precision
    random
    optim
end
    
methods
function obj = MonteCarloPricerSettingsClass(varargin)
%
% 	Constructor. 
% 
% -----------------------------------------------------------------------------------------------
%  Parameters
% -----------------------------------------------------------------------------------------------
%  * Inputs must be given in name-value pairs corresponding to the object properties. If some 
%    properties are not set they may be set to some default values.
%
% -----------------------------------------------------------------------------------------------
%  Output
% -----------------------------------------------------------------------------------------------
%   [1x1 MonteCarloPricerSettingsClass] The object.
% -----------------------------------------------------------------------------------------------

    obj.ParseConstructorInputs(varargin{:});

    % Default values:
    if isempty(obj.tn) && size(obj.n,1)==1;obj.tn=Inf;end
    if isempty(obj.precision);obj.precision = 'double';end

end
function [idxGrp,nActual,nPerYearActual,nPerYearDesired] = GetNumSteps(obj,T)
%
%   Returns the number of simulation steps needed to price each of the input expiries 
%   (and some other information). If Richard-Romberg (RR) extrapolation is enabled we return
%   the number of simulation steps required for both levels of precision as needed by
%   the extrapolation. Outputs 'nActual', 'nPerYearActual' are then [Nx2] sized matrices,
%   and the first columns refer to the finest discretisation level, the second columns to 
%   the coarser discretisation level.
%
% -----------------------------------------------------------------------------------------------
%  Parameter
% -----------------------------------------------------------------------------------------------
%   T: [Nx1 real] Expiries to be priced.
%
% -----------------------------------------------------------------------------------------------
%  Outputs
% -----------------------------------------------------------------------------------------------
%   idxGrp:          [Nx1 integer]          Vector with indices grouping the expiries according 
%                                           to obj.tn.
%
%   nActual:         [Nx1 or Nx2 integer]   The actual/required number of steps to simulate till 
%                                           each expiry.
%
%   nPerYearActual:  [Nx1 or Nx2 real]      The actual/required number of steps per year to 
%                                           simulate till each expiry. May not necessarily be 
%                                           integer valued. Values are always larger than those 
%                                           in nPerYearDesired but usually only differ very 
%                                           slightly (if at all).
%
%   nPerYearDesired: [Nx1 integer]          The desired number of steps per year according to
%                                           obj.tn and obj.n, and that at the finest 
%                                           discretisation level only (if using RR extrapolation).
%
% -----------------------------------------------------------------------------------------------
    
   % Determine which group each expiry belongs to:
   idxGrp = size(obj.tn,1) - sum(bsxfun(@(x,y)(x<=y),T,obj.tn'),2) + 1;
   uniqGrps = unique(idxGrp);
   nGrps = size(uniqGrps,1);   
   
   % Extract the desired number of steps per year (at the finest discretisation level):
   nPerYearDesired = obj.n(idxGrp);

   % Compute total number of steps (and step size) per expiration group:
   [totalStepsPerGrp,dtPerGrp] = deal(NaN(nGrps,1));
   for i=1:nGrps
       idxGrp_i = idxGrp == uniqGrps(i);
       maxT = max(T(idxGrp_i));
       totalStepsPerGrp(i) = floor2(maxT*obj.n(uniqGrps(i)));
       dtPerGrp(i) = maxT/totalStepsPerGrp(i);
   end

   % Compute actual number of steps per expiry:
   [nPerYearActual,nActual] = deal(NaN(size(T,1),1));
   for i=1:nGrps
        idxGrp_i = idxGrp == uniqGrps(i);
        % Expiries that do not fit into the discretisation are approximated by the
        % simulated values on the grid point just below:
        nActual(idxGrp_i,1) = floor2(T(idxGrp_i)./dtPerGrp(i));
        nPerYearActual(idxGrp_i,1) = 1./dtPerGrp(i);
   end
   
end
function ClearNumbers(obj)
%
% 	Clears random numbers.
%
    obj.random = [];
end    
end

end

