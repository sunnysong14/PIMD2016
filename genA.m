function A = genA(nS, seed)
% Usage: 
%   Generate mixing matrix $A$
%%
rng(seed); % contrl random

% Generate A: random with ||column||=1 & plus identity matrix
A = rand(nS,nS);
for ii = 1:nS
    A(:,ii) = A(:,ii)/norm(A(:,ii));
end
A = A + eye(nS);
end  % End OF FUNCTION