function ret = SINR_mean(A, W, Sigma)
% function ret = SINR_mean(A, W, Sigma)
% Inputs:
%   A -- The true mixing matrix from the ICA model.
%   W -- The approximated demixing matrix
%   Sigma -- the Covariance of the Gaussian Noise matrix (before
%   demixing)
% return:
%   SINR_mean -- Returns the mean signal to interference plus noise ratio
%   of the demixed ICA signals.
    
    d = size(A, 1);
    SINR = zeros(d, 1);
    for k = 1:d
        SINR(k) = SINR_Theoretical(A, W, Sigma, k);
    end

    ret = sum(SINR) / d;
end