function [c,alpha,f,K] = ExpandProductOfKernels(K1,K2)
%
% 	Computes weights, fractional indices and non-singular functions for a product of two
%   kernel functions. We then have
%
%       K1(t)*K2(t) = c_1*(t^alpha_1)*f_1(t) + ... + c_n*(t^alpha_n)*f_n(t)
%
%   for some n terms.
%
% -----------------------------------------------------------------------------------------------  
%  Parameters
% -----------------------------------------------------------------------------------------------  
%   K1:   [1x1 KernelFunctionClass] See main description.
%   K2:   [1x1 KernelFunctionClass] See main description.
%
% -----------------------------------------------------------------------------------------------  
% Output
% -----------------------------------------------------------------------------------------------  
%   c:              [nx1 real]                  Weights
%   alpha:          [nx1 real]                  Fractional indices. 
%   f:              [nx1 cell]                  Non-singular functions. Callable as f{i}(t), 
%                                               i = 1,...,n.
%   K:              [1x1 KernelFunctionClass]   K1(t)*K2(t) as a KernelFunctionClass. It is 
%                                               however important to note that the elements of
%                                               K.f may depend on the properties of K1 and K2.
%
% -----------------------------------------------------------------------------------------------      

    % Unpack inputs:
    c1 = K1.c;
    c2 = K2.c;
    alpha1 = K1.alpha;
    alpha2 = K2.alpha;
    f1 = K1.f;
    f2 = K2.f;
    n1 = size(c1,1);
    n2 = size(c2,1);

    % Loop over all combinations of terms from either kernel:
    n = n1*n2;
    [c,alpha] = deal(zeros(n,1));
    [f,f_tilde] = deal(cell(n,1));
    k = 1;
    for i=1:n1
        for j=1:n2
            c(k) = c1(i)*c2(j);
            alpha(k) = alpha1(i) + alpha2(j);

            f1_tmp = @(t)(f1{i}(K1,t));
            f2_tmp = @(t)(f2{j}(K2,t));

            f{k} = @(t)(f1_tmp(t).*f2_tmp(t));
            f_tilde{k} = @(Kobj,t)(f1_tmp(t).*f2_tmp(t));

            k = k + 1; 
        end
    end
    
    % Construct kernel function:
    K = KernelFunctionClass(c,alpha,f_tilde,min(K1.delta,K2.delta));

end
