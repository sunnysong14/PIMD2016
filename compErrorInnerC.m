function error_innerC = compErrorInnerC(A, C)
% Usage:
%   computer the error of innder product matrix C
% 
% verify (Ai, Aj)_C = 0
%   see: matrix DotC and array NormAcol
%   [Criteria] o off(DotC)  should equal zero.
%              o diag(DotC) should equal ||A_i||^2
% 
% Liyan for nips16 on 05-10-2016

%%
m = size(A, 1);
DotC = zeros(m);
for i1 = 1 : m
    Ai = A(:, i1);
    for i2 = i1 : m
        Aj = A(:, i2);
        % verify
        DotC(i1,i2) = Ai'*C*Aj;
    end
end

% measurement: sum(|off(DotC)|), smaller the better.
diagDoc = diag(diag(DotC));

% normalized error:
error_innerC = sum(sum( abs(DotC-diagDoc))) / sum(diag(abs(diagDoc)));
% error_innerC = sum(sum( abs(DotC-diagDoc) ));
end % END OF FUNCTION
