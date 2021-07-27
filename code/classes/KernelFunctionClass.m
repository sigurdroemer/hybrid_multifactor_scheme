classdef KernelFunctionClass < GenericSuperClass & dynamicprops
%
%   Holds a kernel function K(t) of the form
%
%       K(t) = c_1*(t^alpha_1)*f_1(t) + ... + c_n*(t^alpha_n)*f_n(t)                                                            
%
%   where (c_1,...,c_n) are non-negative numbers, (alpha_1,...,alpha_n) fractional indices, and 
%   (f_1,...,f_n) are non-singular (non-negative) functions. It is assumed (c,alpha,f) are chosen
%   so f is completely monotone.
%
% -----------------------------------------------------------------------------------------------
%  Properties  
% -----------------------------------------------------------------------------------------------
%   * Remark: More parameters than the below can be added. See inputs 'addParNames' and 
%     'addParVals' under the class constructor.
%
%   c:              [nx1 real]      Weights
%   alpha:          [nx1 real]      Fractional indices. 
%   f:              [nx1 cell]      Non-singular functions. Each element is of the form 
%                                   [1x1 function] and must be callable as either f(obj,t), t >= 0, 
%                                   where obj refers to the particular KernelFunctionClass object. 
%                                   Must also be vectorised in t.
%
%   delta:          [1x1 real]      A very small positive number used for computing various types 
%                                   of integrals involving K. Used by methods Integrate and
%                                   IntegrateVolterra. Is then passed directly to 
%                                   function IntegrateFractional when used. 
%
% -----------------------------------------------------------------------------------------------

properties
    c
    alpha
    f
    delta
end

methods
function obj = KernelFunctionClass(c,alpha,f,delta,addParNames,addParVals)
%
%   Constructor.
%
% -----------------------------------------------------------------------------------------------
%  Parameters
% -----------------------------------------------------------------------------------------------  
%  * See class description for an explanation of c, alpha, f, delta. Remaining inputs are 
%    explained below (can be left empty).
%
%   addParNames:   [1xN cell] Each element is a string representing a variable name that will 
%                             then be created in the object. Intended use: Allows obj.f to depend 
%                             on other parameters besides obj.alpha.
%
%   addParVals:    [1xN cell] Values to initialise each additional parameter at.
%
% -----------------------------------------------------------------------------------------------  

    if exist('addParNames','var')
        for i=1:size(addParVals,2)
            obj.addprop(addParNames{i});
            obj.(addParNames{i}) = addParVals{i};
        end
    end

    if exist('c','var')  && ~isempty(c)
        obj.c = c;
        if size(c,2) > 1
            error('KernelFunctionClass: ''c'' must be a column vector.');
        end
    end    
    
    if exist('alpha','var')  && ~isempty(alpha)
        obj.alpha = alpha;
        if size(alpha,2) > 1
            error('KernelFunctionClass: ''alpha'' must be a column vector.');
        end        
    end

    if exist('f','var') && ~isempty(f)
        if size(f,2) > 1
            error('KernelFunctionClass: ''f'' must be a column cell array.');
        end      
        if ~iscell(f) && size(f,1) > 1
            error('KernelFunctionClass: ''f'' must be a cell array.');
        elseif ~iscell(f) && size(f,1) == 1
            obj.f = {f};
        else
            obj.f = f;
        end
    end
    
    if exist('delta','var') && ~isempty(delta)
        if delta < 0
            error('KernelFunctionClass: ''delta'' must be non-negative.');
        end     
        obj.delta = delta;
    end        

end
function val = Eval(obj,t)
%
% 	Evaluate the kernel.
%
% -----------------------------------------------------------------------------------------------  
%  Parameters
% -----------------------------------------------------------------------------------------------  
%   t:   [NxM real] Input values.
%
% -----------------------------------------------------------------------------------------------  
% Output
% -----------------------------------------------------------------------------------------------  
%   val: [NxM real] Output values.
%
% -----------------------------------------------------------------------------------------------  
    
    val = 0;
    for i=1:size(obj.c,1)
        val = val + obj.c(i)*(t.^obj.alpha(i)).*obj.f{i}(obj,t);
    end

end
function val = Integrate(obj,a,b,p,nu,c,gamm)
%
%   Uses numerical integration to evaluate an integral of the form 
%
%       int_a^b K(s)^p ds                                                        
%
%   where p = 1 or 2 and 0 <= a <= b.
%
% -----------------------------------------------------------------------------------------------  
%  Parameters
% -----------------------------------------------------------------------------------------------  
%   a:      [Nx1 real]      See the description.
%   b:      [Nx1 real]      See the description.
%   p:      [1x1 real]      See the description.
%
%   The following additional parameters also apply:
%
%   nu:     [1x1 real]      To compute the integral we approximate K(s)^p, s > 0, exactly for
%                           0 < s <= nu and by a sum of exponentials as
%
%                               K(s)^p approx(=) sum_{j=1}^m c_j exp(-gamm_j*s)
%
%                           for s > nu. 
%
%                           IMPORTANT: Note that the sum of exponentials approximation must be 
%                           for K "raised to the power p" and not for K itself. 
%                           
%                           Note also:
%                           The sum of exponentials approximation is justified since the product 
%                           of completely monotone functions is also completely monotone.
%
%   c:      [1xm real]      Coefficients for sum of exponentials approximation.
%   gamm:   [1xm real]      Exponents for sum of exponentials approximation.
%
% -----------------------------------------------------------------------------------------------  
%  Output
% -----------------------------------------------------------------------------------------------  
%   val:    [Nx1 real] Value of integral(s).
%
% -----------------------------------------------------------------------------------------------      

    if size(a,1) ~= size(b,1) || size(a,2) > 1 || size(b,2) > 1
        error('KernelFunctionClass:Integrate: ''a'' and ''b'' must be of size Nx1.');
    end
    
    if any(b < a | a < 0)
        error('KernelFunctionClass:Integrate: We require 0 <= a <= b.');
    end
    
    if size(c,1) > 1 || size(gamm,1) > 1 
        error('KernelFunctionClass:Integrate: ''c'' and ''gamm'' must be row vectors.');
    end    

    [I1,I2] = deal(zeros(size(b)));
    idxBelowNu = a < nu;
    idxAboveNu = b > nu;

    %% Compute the singular part:
    % Remark: A representation similar to 
    %   K(t) = c_1*(t^alpha_1)*f_1(t) + ... + c_n*(t^alpha_n)*f_n(t),
    % holds for K(t)^2 as well.
    
    if p == 1
        K_dummy = obj;
    elseif p == 2
        [~,~,~,K_dummy] = ExpandProductOfKernels(obj,obj);
    end    
    c_weight = K_dummy.c;
    alpha = K_dummy.alpha;
    f = K_dummy.f;
    
    idx = find(idxBelowNu);
    for i=1:size(idx,1)
        for j=1:size(c_weight,1)
            I1(idx(i)) = I1(idx(i)) + c_weight(j)*IntegrateFractional(alpha(j),@(s)(f{j}(obj,s)),...
                                                            a(idx(i)),min(nu,b(idx(i))),obj.delta);
        end
    end
    
    %% Compute the non-singular part:
    if any(idxAboveNu)
        idxZero = gamm == 0;
        if any(~idxZero)
            I2(idxAboveNu) = sum((c(~idxZero)./gamm(~idxZero))...
                                 .*(exp(-max(nu,a(idxAboveNu))*gamm(~idxZero)) ...
                                 - exp(-b(idxAboveNu)*gamm(~idxZero))),2);                 
        end
        if any(idxZero)
            I2(idxAboveNu) = sum(c(idxZero)).*(b(idxAboveNu)-max(nu,a(idxAboveNu)));
        end
    end
    
    % Combine results:
    val = I1 + I2;


end
function val = IntegrateVolterra(obj,a,b,c,p,g)
%
%   Uses numerical integration to evaluate an integral of the form 
%
%       int_a^b K(c-s)^p g(s) ds                                                        (*)
%
%   where p = 1 or 2, 0 <= a <= b <= c and g is a non-singular function on [a,b]. 
%
%   A note on the implementation: We start by defining d = max(c-obj.delta,a). If d >= b we then 
%   approximate (*) using Matlab's 'integral' function. If d < b we instead decompose the integral 
%   as
%
%       int_a^b K(c-s)^p g(s) ds = I_1 + I_2
%
%   where
%
%       I_1 := int_a^d K(c-s)^p g(s) ds
%   
%       I_2 := int_d^b K(c-s)^p g(s) ds.
%
%   I_1 is then computed using Matlab's 'integral' function while I_2 is computed using the 
%   'IntegrateFractional' function.
%
% -----------------------------------------------------------------------------------------------  
%  Parameters
% -----------------------------------------------------------------------------------------------  
%   a:      [1x1 real]      See the description.
%   b:      [1x1 real]      See the description.
%   c:      [1x1 real]      See the description.
%   p:      [1x1 real]      See the description.
%   g:      [1x1 function]  See the description. Must be vectorised.
%
% -----------------------------------------------------------------------------------------------  
%  Output
% -----------------------------------------------------------------------------------------------  
%   val:    [1x1 real]      Value of integral.
%
% -----------------------------------------------------------------------------------------------  
    
    if size(a,1) > 1 || size(b,1) > 1 || size(a,2) > 1 || size(b,2) > 1
        error('KernelFunctionClass:Integrate: ''a'' and ''b'' must be of size 1x1.');
    end
    
    if a == b
        val = 0;
        return;
    end
    
    if ~(0 <= a && a <= b && b <= c)
        error('KernelFunctionClass:Integrate: We require 0 <= a < b <= c.');
    end    

    d = max(c-obj.delta,a);
    
    if d >= b
        val = integral(@(s)((obj.Eval(c-s).^p).*g(s)),a,b);
        return;
    end
    
    % Here: d < b:
    
    %% Compute I_1:
    I_1 = integral(@(s)((obj.Eval(c-s).^p).*g(s)),a,d);

    %% Compute I_2. 
    % Note that
    %   I_2 = int_d^b K(c-s)^p g(s) ds = int_{c-b}^{c-d} K(s)^p g(c-s) ds,
    %
    % and that a representation similar to 
    %
    %   K(t) = c_1*(t^alpha_1)*f_1(t) + ... + c_n*(t^alpha_n)*f_n(t),
    %
    % holds for K(t)^2.
    
    if p == 1
        K_dummy = obj;
    elseif p == 2
        [~,~,~,K_dummy] = ExpandProductOfKernels(obj,obj);
    end    
    c_weights = K_dummy.c;
    alpha = K_dummy.alpha;
    f = K_dummy.f;
    
    I_2 = 0;
    for i=1:size(c_weights,1)
        I_2 = I_2 + c_weights(i)*IntegrateFractional(alpha(i),@(s)(f{i}(obj,s).*g(c-s)),...
                                                     c-b,c-d,obj.delta);
    end
    
    val = I_1 + I_2;

end
function plot(obj,Ks,a,b,nPts)
%
%   Creates a simple plot showing the function.
%
% ---------------------------------------------------------------------------------------------
% Parameters 
% ---------------------------------------------------------------------------------------------
%   Ks:      [1xN cell]      Each element is a [1x1 KernelFunctionClass] that will then be 
%                            plotted. If empty we default to only plotting the current object.
%   a:       [1x1 real]      Lower end-point. Default is 1/1000.
%   b:       [1x1 real]      Upper end-point. Default is 1.
%   nPts:    [1x1 integer]   Number of points to sample between a and b. Default is 1000.
%
% ---------------------------------------------------------------------------------------------
%    
    
    if ~exist('a','var') || isempty(a)
        a = 1/1000;
    end
    if ~exist('b','var') || isempty(b)
        b = 1;
    end
    if ~exist('nPts','var') || isempty(nPts)
        nPts = 1000;
    end    
    
    dx = (b-a)/(nPts-1);
    x = (a:dx:b)';
    
    if ~exist('Ks','var') || isempty(Ks)
        Ks = {obj};
    end
    
    figure;
    for i=1:size(Ks,2)
        y = Ks{i}.Eval(x);
        plot(x,y,'DisplayName',['Curve ',num2str(i)]);hold on;
    end
    xlabel('t');
    ylabel('Values');
    legend();
    xlim([a-0.1,b+0.1]);
    hold off;
    
end
end

end

