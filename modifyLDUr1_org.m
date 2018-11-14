function [Lm, Dm, Um] = modifyLDUr1(L, D, U, Ml, a)
% [Lm, Dm, Um] = calcLDU(L*D*U + Ml*a*Ml');
N = size(L, 1);
for i = 1:N-1
    alpha = a * Ml(i);
    D(i, i) = D(i, i) + alpha * Ml(i);
    beta = alpha / D(i, i);
    Ml(i+1:end) = Ml(i+1:end) - U(i, i+1:end).' * Ml(i);
    U(i, i+1:end) = U(i, i+1:end) + beta * Ml(i+1:end).';
    a = a - alpha * beta;
end
alpha = a * Ml(N);
D(N, N) = D(N, N) + alpha * Ml(N);
Lm = U.';
Um = U;
Dm = D;
end