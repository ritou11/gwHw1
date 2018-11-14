function [Lm, Dm, Um] = modifyLDUlr(D, U, dA)
% [Lm, Dm, Um] = calcLDU(U.'*D*U + dA);
N = size(U, 1);
Sm = generateSm(unique(floor((find(dA) - 1) / N + 1)), U);
U22 = U(Sm, Sm);
D22 = D(Sm, Sm);
A22p = U22.' * D22 * U22;
A22ps = A22p + dA(Sm, Sm);
[~, D22s, U22s] = calcLDU(A22ps);
D(Sm, Sm) = D22s;
U(Sm, Sm) = U22s;
Lm = U.';
Dm = D;
Um = U;
end