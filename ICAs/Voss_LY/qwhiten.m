% Liyan adjusts the format and add comments simply. Not a big deal. 
% Records of changing:
%   05-08-2016: codes on- & off- centering
%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Y, W, b] = qwhiten(X, op)
%   function [Y, W, b] = qwhiten(X, op)
%   qwhiten is short for quasi-whitening.
%
%   Inputs:
%   X -- The inputted data.  Each column of the matrix is a random
%   sample of the linearly mixed observations.  The number of columns gives
%   the number of samples.
%
%   op -- 'whiten', 'id quas-orth'.  The "identity" based
%   algorithm is deterministic in nature, but requires more operations
%   than "random" (the time difference is a constant factor).  Both
%   "identity" and "random" are quasi-whitening algorithms robust to
%   Gaussian noise, whereas "whiten" is a standard whitening algorithm
%   based on PCA that is not robust to Gaussian noise.
%   "half rnd" only works when the signals all have kurtosis of the same
%   sign.
%
%   returns:
%       Y -- The quasi-orthogonalized (or whitened) data.
%       W -- The resulting quasi-whitening (or whitening) matrix.
%       b -- The offset of the data from the origin before preprocessing.
%   The samples in Y are computed as y = W(x-b)

%%
[d, n] = size(X);    
C = zeros(d);

%% LY TMP 05-06-2016
% Note: TAG_CENTER is a global variable, that are only used in pegi
% and our method. 

global TAG_CENTER

% LY centering
if TAG_CENTER == 1
    X_orig = X;  % LY: keep original X for analysis
    b = mean(X, 2);
    X = bsxfun(@minus, X, b);  % X = X - repmat(b, 1, N);
else
    % LY, TMP add for correct return
    b = zeros(size(X, 1), 1);
end

%%
if strcmp( op, 'whiten' )
    C = cov(X');
    W = C^(-0.5);
    Y = W * X;
    return
end

% Generate the initial hessian / cumulant tensor matrix.
XXt = X*X';  % Compute in order to (unnoticeably) increase hessian efficiency.
if strcmp( op, 'id quas-orth' )
    for i = 1:d
        u = canonvec(i, d);
        C = C + cum4hes(X, u, XXt);
    end
    C = C / 12;  % Division by 12 makes this equivalent to a change of variable in the fourth cumulant tensor techniques.
else
    fprintf(2, ['ERROR:  Invalid option flag:  ' op]);
    return
end

% Inversion and second step of the quasi-whitening.
C2 = zeros(d);
[U, D] = eig(C);
% will be using the inverse of C, and placing half of the effect of D
% in the U vectors.

for i = 1:d
    C2 = C2 + (1/D(i, i))*cum4hes(X, U(:, i), XXt);
end
C2 = C2 / 12;
[U, Lambda] = eig(C2);
len = size(Lambda, 1);
Lambda_mask = (diag(Lambda) > 0);
if (sum(Lambda_mask) < len)
   warning(['Decomposition matrix for quasi-orthogonalization is not positive definite.  ' ...
       'Ignoring negative eigenvalue directions.  Output could very easily be bogus.  ' ...
       'Standard whitening may provide better results.']);
   Lambda;
   Lambda(1:len+1:end) = Lambda * Lambda_mask;
   W = pinv(Lambda)^(1/2) * U';
else
   W = Lambda^(-1/2)*U';
end
Y = W*X;
end  %% END OF FUNCTION

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sub-functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ei = canonvec(i, dim)
%   function ei = canonvec(i, dim)
%   Produces the ith canonical vector in a dim dimensional space.
    ei = zeros(dim, 1);
    ei(i) = 1;
end

function u = rndunitvec(dim)
%  Generates a random vector in dim uniformly from the unit sphere.
    u = normrnd( zeros(dim, 1), ones(dim, 1) );
    u = u / norm(u);
end
