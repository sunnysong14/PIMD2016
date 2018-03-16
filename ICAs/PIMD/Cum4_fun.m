function Cum4_tsor = Cum4_fun(X)
% Usage:
%   Compute the sample the 4th-order cumulants. 
%   The codes were originally from JADE method that were commented off to improve the efficiency.
% Inputs:
%   X -- centered X with the size d x N, where d is the dimension and N is
%   the samplings.
%   The row corresponds each source, and the columns are the samplings of each source.
% Output:
%   Cum4_tsor -- Q_x(i,j,k,l), the 4th-order cumulants with the dimension n^4.
% 
% [Note] TO-DO
%   Low efficiency. Pls upspeed it later. Refer JADE's codes.
% 
% Liyan for NIPS16 05-04-2016

%%
m = size(X, 1);  % #row should equal to the source dimension
if m >= size(X,2)
    warning('Warning: note the form of observation matrix X: the rows should corresponding to each source, and columns to samplings.');
end

% should I do centering as kurtosis()? The codes below are from kurtosis()
% ----------------------------------------
% % Center X, compute its fourth and second moments, and compute the
% % uncorrected kurtosis.
% x0 = x - repmat(nanmean(x,dim), tile);
% s2 = nanmean(x0.^2,dim); % this is the biased variance estimator
% m4 = nanmean(x0.^4,dim);
% k = m4 ./ s2.^2;

% main calculation
Cum4_tsor = zeros(m, m, m, m) ;  % store the 4th order cumulants
for i1 = 1 : m
    for i2 = 1 : m
        for i3 = 1 : m
            for i4 = 1 : m
                % the 4th-order cumulant formula
                Cum4_tsor(i1,i2,i3,i4) = ...
                    mean(X(i1,:).* X(i2,:).* X(i3,:).* X(i4,:)) ...
                    - mean(X(i1,:).* X(i2,:)) * mean(X(i3,:).* X(i4,:))...
                    - mean(X(i1,:).* X(i3,:)) * mean(X(i2,:).* X(i4,:)) ...
                    - mean(X(i1,:).* X(i4,:)) * mean(X(i2,:).* X(i3,:));
            end
       end
    end
end
end%END OF FUNCTION
