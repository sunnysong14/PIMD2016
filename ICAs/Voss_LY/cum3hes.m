% function H = cum3hes(X, u)
% H -- Returned hessian matrix.
% X -- Matrix where each column represents a data point.  The data is
%      assumed to be centered (i.e., mean subtracted).
% u -- directional vector.  The Hessian of the u-directional fourth cumulant
%      is taken.
% Computes and returns an approximation to the hessian with respect to the
% directional variable u of k_3(u'x) where x denotes the random variable
% being sampled.  Here, k_3 is the fourth k-statistic estimate to the
% fourth cumulant.

function H = cum3hes(X, u)
    %% compute the unbiased estimator
    [d, n] = size(X);
    
    D = sparse(1:n, 1:n, u'*X);  % generates a sparse diagonal matrix with 
                                 % entries from u' * X
    H = (X * D * X') * (6*n)/((n-1)*(n-2));
end
