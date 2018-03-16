function QxM_mat = QxM_fun(Cum4_tsor, M)
% Usage:
%   Compute the equation $Q_x(M)$
%   The codes were originally from JADE method that were commented off to improve the efficiency.
% 
% Input:
%   Cum4_tsor -- the 4th-order cumulant tensor
%   M -- the arbitrary matrix on which the 4th-order cumulant tensor is
%   operated on.
% Output:
%   QxM_mat -- the matrix of $Q_x(M)$
% 
% Liyan for NIPS16 05-04-2016

%%
m = size(Cum4_tsor, 1);
if m ~= length(M)
    error('Error: the dim of the arbitrary matrix M should equal to that of Cum4_tsor.')
end

% main calculation
QxM_mat = zeros(m, m);
for i1 = 1 : m
    for i2 = 1 : m
      QxM_mat(i1, i2) = sum(sum(squeeze(Cum4_tsor(i1,i2,:,:)) .* M));
    end
end
end  % END OF FUNCTION
%% test
%{
m = size(Cum4_tsor, 1);
M = rand(m);

% main call
QxM_mat = QxM_fun(Cum4_tsor, M);
%}
