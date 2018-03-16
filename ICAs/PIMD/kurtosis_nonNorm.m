function nfKurt_arr = kurtosis_nonNorm(X)
% Usage:
%   Compute non-uniformed kurtosis with Eq.(2.90) [Aapo_ICA2001,P38].
%   Thus, the minimum kurtosis is not -2.
%   The reason is to verify the following equation: for any matrix M
%   $Q_x(M) = A*D*A^T$, where A is the mixing matrix, and D is diagonal
%   with $D_{ii} = k_4(s_i)*A_i^T*M*A_i$, A_i is the i-th column of A.
% Input:
%   X -- each row corresponds with one dimension, and each column contains
%   one sampling.
% 
% Liyan for NIPS16, 05-04-2016

%%
m = size(X,1);  % dimension
nfKurt_arr = -ones(m,1);

for d = 1 : m
    Xi = X(d,:);
    nfKurt_arr(d) = mean(Xi.^4)- 3*mean(Xi.^2)^2;  % non-normalized form
end

end  % END OF FUNCTION
