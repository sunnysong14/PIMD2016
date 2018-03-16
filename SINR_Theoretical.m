% Liyan integrates this function from Voss nips15 codes.
% 04-20-2016 for nips16.
%%
function SINR = SINR_Theoretical(A, W, Sigma, k)
% function SINR = SINR_Theoretical(A, W, Sigma, k)
% Inputs:
%   A -- The true mixing matrix from the ICA model.
%   W -- The approximated demixing matrix
%   Sigma -- the Covariance of the Gaussian Noise matrix (before
%   demixing)
%   k -- Signal of which to take the SINR.  k should be at most the
%   ambient dimension.
% return:
%   SINR_Theoretical -- Rference plus noise ratio
%   of the demixed signals.eturn the kth signal to inter
%   SINR is given from the actual (not sample) variance of the signal
%   divided by the actual variance of the noise plus interference.
%   Further, the convention is used that if 0 signal variance is recovered,
%   then the resulting signal to interference ratio is 0 regardless of whether
%   the denominator is 0.
    
    %% Calculate the approximately demixed signal contributions.
    b = W(k, :) * A;
    
    %% Relevant Variances
    % [sig_var, index] = max(b.^2);  % We assume that we recovered the signal that has the maximum variance.
    sig_var = b(k)^2;
    if sig_var == 0
        SINR = 0;
    else
        interference = norm(b)^2 - sig_var;
        noise = W(k, :) * Sigma * W(k, :)';

        %% SINR calculation
        SINR = sig_var / ( interference + noise );
    end
end