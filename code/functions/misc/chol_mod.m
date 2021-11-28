function A = chol_mod(SIGMA,tol)
% 
%   Performs Cholesky factorization of a covariance matrix. Even works for positive semi-definite 
%   matrices in which case we use the Modified Cholesky Algorithm of (Cheng and Higham, 1998) to
%   slightly pertube the input covariance matrix.
%
% -------------------------------------------------------------------------------------------------
%  Parameters
% -------------------------------------------------------------------------------------------------
%   SIGMA:   [NxN real]           Covariance matrix; must be positive semi-definite.
%   tol:     [1xM real or empty]  Tolerance for checking eigenvalues is norm(SIGMA,'fro')*sqrt(tol).
%                                 If 'tol' is a vector, code is run from the first to the last 
%                                 entry until the factorization succeeds. Default is a [1x1 real] 
%                                 variable equal to eps.
%
% -------------------------------------------------------------------------------------------------
%  Output
% -------------------------------------------------------------------------------------------------
%   A:       [NxN real]           Lower diagonal matrix such that A*transpose(A) = SIGMA.
%
% -------------------------------------------------------------------------------------------------
%  References
% -------------------------------------------------------------------------------------------------
%   o Cheng, S.H., and Higham, N.J., A modified Cholesky algorithm based on a symmetric indefinite 
%     factorization. SIAM Journal on Matrix Analysis and Applications, 1998, 19(4), 1097-1110.
%
% -------------------------------------------------------------------------------------------------

    % Check if matrix is symmetric:
    if ~issymmetric(SIGMA)
        error('chol_mod: Covariance matrix must be symmetric.');
    end

    % Check if matrix is positive semidefinite within tolerance:
    if ~exist('tol','var') || isempty(tol);tol=eps;end
    
    for i=1:size(tol,2)
        try
            tol_eig = norm(SIGMA,'fro')*sqrt(tol(i));
            if any(eig(SIGMA) < -tol_eig)
                error('chol_mod: Matrix was not deemed positive semidefinite within tolerance.');
            end

            try
                % Attempt a Cholesky factorisation:
                A = chol(SIGMA, 'lower');
                return;
            catch
                % If the matrix is not positive definite:
                [L, D, P] = modchol_ldlt(SIGMA,tol_eig);
                SIGMA_pert = P'*L*D*L'*P;
                A = chol(SIGMA_pert, 'lower');
                return;
            end
        
        catch
            % Try with a new tolerance value.
            
        end
    
    end
    
    error('chol_mod: Matrix was not deemed positive semidefinite within tolerance(s).');
    
end

