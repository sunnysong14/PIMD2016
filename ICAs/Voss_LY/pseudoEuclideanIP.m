function M = pseudoEuclideanIP(X, cumulant)
%   function M = pseudoEuclideanIP(X, op)
%   Inputs:
%   X -- The inputted data matrix.  Each column of the matrix is a random
%        observation of the ICA model.  The number columns gives the number
%        of samples.
%   cumulant -- One of the following:
%        'k3' -- Third cumulant based construction of the pseudo-Euclidean
%            inner product space.
%        'k4' -- Fourth cumulant based construction of the pseudo-Euclidean
%            inner product space.
%   Returns:
%   M -- A matrix approximately of the form inv(ADA') where D is a diagonal
%        matrix with possibly negative entries.  In particular, the matrix
%        M forms a psuedo-inner product for a space in which the columns of
%        A are mutually orthogonal.
    
    [d, n] = size(X);
   
    %% Determine which cumulant to use.
    if strcmp(cumulant, 'k3')
        cumhes = @( u ) (1/6) * cum3hes(X, u); % Division by 1/6 makes this mathematically equivalent to methods based on collapsing a 3rd k-statistic tensor.
    elseif strcmp(cumulant, 'k4')
        XXt = X*X'; % Computed for efficiency reasons in the cumulant 
                    % hessian operation.
        cumhes = @( u ) (1/12) * cum4hes(X, u, XXt); % Division by 1/12 makes this mathematically equivalent to methods based on collapsing the 4th k-statistic tensor.
    else
        error('Invalid choice of cumulant operation:  ', cumulant);
    end
    
    %% Construct the inner product
    M = zeros(d, d);
    for i = 1:d;
        % construct the ith canonical vector
        u = sparse([i], [1], [1], d, 1);  % LY: e_i simply
        
        % compute the corresponding hessian, and accumlate within the
        % desired (not yet inverted) inner product matrix contruction.
        M = M + cumhes(u);
    end
    M = inv(M);
end
