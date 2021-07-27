function W_bold = SampleVolterra(A,K,kappa,dt,antithetic,N,n,Y)
%
%   Sample (correlated) Gaussian Volterra integrals. 
%
%   Let A be a mxm matrix such that rho := AA' is a valid correlation matrix. Let then 
%   Z(t) = (Z_1(t),...,Z_m(t))' be a independent Brownian vector and define W(t) = A*Z(t).
%   For i=1,...,m, let kappa_i be non-negative integer and K_i be completely monotone kernel 
%   function. 
%
%   The code then samples independent realisations of the following vector:
%
%                W_1(dt) - W_1(0), 
%                int_{0}^{dt} K_1(dt-s) dW_1(s),
%                int_{0}^{dt} K_1(2*dt-s) dW_1(s),
%                ...,                                                  
%                int_{0}^{dt} K_1(kappa_1*dt-s) dW_1(s)
%
%                ...,
%
%                W_m(dt) - W_m(0), 
%                int_{0}^{dt} K_m(dt-s) dW_m(s),
%                int_{0}^{dt} K_m(2*dt-s) dW_m(s),
%                ...,                                                  
%                int_{0}^{dt} K_m(kappa_m*dt-s) dW_m(s).
%
%   Remark: Intended use is for a simulation algorithm; this is reflected in the way the function 
%   is set up and the terminology that is used.
%
% ------------------------------------------------------------------------------------------------
%  Parameters
% ------------------------------------------------------------------------------------------------
%   A:            [mxm real]        Factorization of correlation matrix rho.
%   K:            [1xm cell]        Kernel functions. Each element must be an object of class
%                                   [1x1 KernelFunctionClass].
%   kappa:        [1xm integer]     See the main description.
%   dt:           [1x1 real]        See the main description.
%   antithetic:   [1x1 logical]     Also generates antithetic samples (we return realizations 
%                                   for -W(t) too).
%   N:            [1x1 integer]     Number of paths. Must be divisible by 2 if antithetic = true.
%   n:            [1x1 integer]     Number of steps.
%   Y:            [1xm cell]        Underlying N(0,1)'s. The i'th element is used to generate
%                                   terms related to Z_i(t). Can be left empty in which case 
%                                   the random variables are automatically sampled in the code. 
%                                   If non-empty the following applies:
%
%                                   The i'th element should take the form [d(i) x Nindep*n real] 
%                                   where d(i) >= sum(kappa(~idxZero)) + 1 with 
%                                   idxZero = A(:,i) == 0. If all(A(:,i)==0) it is however also
%                                   valid to leave it empty.
%
%                                   Nindep = N if antithetic = false, otherwise Nindep = N/2.
%
% ------------------------------------------------------------------------------------------------
%  Output
% ------------------------------------------------------------------------------------------------
%   W_bold: [1xm cell]  The i'th element is of size [N x n x kappa(i) + 1] containing (in the 
%                       order listed in the main description) the stochastic Volterra integrals 
%                       related to the i'th Brownian motion W_i(t). The bottom half of the paths 
%                       will be antithetic if so chosen.
%
% ------------------------------------------------------------------------------------------------

   % Validate inputs:
   m = size(A,1);
   if size(A,2) ~= m || size(K,1) > 1 || size(K,2) ~= m || (exist('Y','var') && ~ isempty(Y) && ...
           (size(Y,1) > 1 || size(Y,2) ~= m))
       error('SampleCorrelatedVolterraIntegrals: Input dimensions are invalid.');
   end
   
   rho = A*A';
   if any(abs(diag(rho)-1) > 10^(-4))
       error(['SampleCorrelatedVolterraIntegrals: AA'' was not deemed a valid ',....
              'correlation matrix within tolerance.']);
   end

   % Check if antithetics are requested:
   if antithetic
       if mod(N,2) ~= 0
           error(['SampleCorrelatedVolterraIntegrals: When ''antithetic'' = true we require ',...
                  'N divisible by 2.']);
       end
       Nindep = N/2;
   else
       Nindep = N;
   end
   
   % Initialize outputs:
   W_bold = cell(m,1);
   for i=1:m
       W_bold{i} = zeros(N,n,kappa(i)+1);
   end
   
   for j=1:m
       % We sample, compute and add terms related to the j'th independent Brownian motion Z_j(t).
       %
       % More precisely, we first sample the following terms
       %
       % Z_j(dt) - Z_j(0)
       %
       % int_{0}^{dt} K_1(dt-s) dZ_j(s),
       % int_{0}^{dt} K_1(2*dt-s) dZ_j(s),
       % ...,                                                  
       % int_{0}^{dt} K_1(kappa_1*dt-s) dZ_j(s)       
       %
       % ...,   
       %
       % int_{0}^{dt} K_m(dt-s) dZ_j(s),
       % int_{0}^{dt} K_m(2*dt-s) dZ_j(s),
       % ...,                                                  
       % int_{0}^{dt} K_m(kappa_m*dt-s) dZ_j(s).
       %
       % We then (appropriately) multiply the j'th column of A onto these samples and add them
       % to the appropriate outputs.
       %
       % To avoid unnecessary computational costs we don't compute terms where the corresponding 
       % entries in A are zero.
       %
       
        % Find the non-zero elements of the j'th column of A:
        idx = A(:,j) ~= 0;
        
        if any(idx)
            % Create and factorize covariance matrix for K's convolved with Z_j(t):
            SIGMA = GetVolterraCovarianceMatrix(K(idx),kappa(idx),dt);
            B = chol_mod(SIGMA,10.^(-(16:-1:10)));

            % Verify input size of Y:
            if exist('Y','var') && ~isempty(Y) && size(Y{j},2) < sum(kappa(idx)) + 1
               error(['SampleCorrelatedVolterraIntegrals: Invalid size of ''Y{',...
                      num2str(j),'}''.']);                
            end
            
            % Sample Volterra terms:
            if ~exist('Y','var') || isempty(Y)
                K_dZj = randn(Nindep*n,sum(kappa(idx))+1)*B';
            else
                K_dZj = Y{j}(1:Nindep*n,1:(sum(kappa(idx))+1))*B';
            end

            % Add Volterra samples to appropriate outputs (after multiplying by A):
            idx = 2;
            for i=1:m
                if A(i,j) ~= 0
                    % Add terms:
                    W_bold{i}(1:Nindep,:,:) = W_bold{i}(1:Nindep,:,:) ...
                                              + A(i,j)*reshape(K_dZj(:,[1,idx:idx+kappa(i)-1]),...
                                                                    Nindep,n,kappa(i)+1);
                    % Update counter:
                    idx = idx + kappa(i);                                                                
                end
            end
        end
   end
   
   % Compute antithetic samples if so desired:
   if antithetic
       for i=1:m
            W_bold{i}(Nindep+1:end,:,:) = -W_bold{i}(1:Nindep,:,:);
       end
   end

end

