function covM = CovMatrixHybrid(n,kappa,alpha,prec)
%
%   Returns the covariance matrix needed for the hybrid scheme of (Bennedsen et al., 2017).
%
%   Remarks:
%       o Function uses Matlab's Symbolic Toolbox.
%
% -------------------------------------------------------------------------------------------------
%  Parameters
% -------------------------------------------------------------------------------------------------
%   n:      [1x1 integer] Number of steps per year.
%   kappa:  [1x1 integer] Number of sub-integrals to be simulated exactly.
%   alpha:  [1x1 real] Roughness index, equals Hurst exponent minus 1/2.
%   prec:   [1x1 string] Precision, options are 'double' and 'single'.
%
% -------------------------------------------------------------------------------------------------
%  Output
% -------------------------------------------------------------------------------------------------
%   covM:   [(kappa+1)x(kappa+1) real] Covariance matrix.
%
% -------------------------------------------------------------------------------------------------
%  References
% -------------------------------------------------------------------------------------------------
%   o Bennedsen, M., Lunde, A. and Pakkanen, M.S., Hybrid scheme for Brownian semistationary 
%     procesess. Finance and Stochastics, 2017, 21(4), 931-965.
%
% -------------------------------------------------------------------------------------------------
%

    covM = nan(kappa + 1, kappa + 1, prec);
    covM(1, 1) = 1 / n;

    for j=2:size(covM, 1)
       covM(1, j) = ( (j-1)^(alpha + 1) - (j-2)^(alpha + 1) ) ...
                    / ((alpha + 1)*n^(alpha + 1));
       covM(j ,j) = ( (j-1)^(2*alpha + 1) - (j-2)^(2*alpha + 1) )  ...
                    / ((2*alpha + 1)*n^(2*alpha + 1));
    end

    for j=2:size(covM, 1)
        for k=2:size(covM, 2)
            if j < k
                covM(j, k) = (1 / ((alpha + 1)*n^(2*alpha + 1)) )*...
                    ( ((j-1)^(alpha + 1))*((k-1)^(alpha))...
                    *hypergeom([-alpha, 1], alpha + 2, (j - 1)/(k - 1)) ...
                    - ((j - 2)^(alpha + 1))*((k-2)^(alpha))...
                    *hypergeom([-alpha, 1], alpha + 2, (j - 2)/(k - 2))  );
            end
        end
    end
    covM(isnan(covM))=0;
    covMDiagZero =  covM - diag(diag(covM));
    covM = covM + triu(covMDiagZero)';
    
end
