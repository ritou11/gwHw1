function [L, D, U] = calcLDU(A)
[L, U] = calcLU(A);
D = diag(diag(L));
N = size(A, 1);
for i=1:N
    L(:, i) = L(:, i) / D(i, i);
end
end