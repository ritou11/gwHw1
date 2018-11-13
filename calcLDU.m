function [L, D, U] = calcLDU(A)
[L, U] = calcLU(A);
D = diag(diag(L));
L(logical(speye(size(L)))) = 1;
end