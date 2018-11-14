function [Lm, Dm, Um] = modifyLDUr1(D, U, Ml, a)
% [Lm, Dm, Um] = calcLDU(L*D*U + Ml*a*Ml');
Sm = generateSm(find(Ml), U);
for idx = 1:size(Sm,1)
    i = Sm(idx);
    alpha = a * Ml(i);
    D(i, i) = D(i, i) + alpha * Ml(i);
    beta = alpha / D(i, i);
    for jdx = 1:size(Sm,1)
        j = Sm(jdx);
        if j <= i
            continue;
        end
        Ml(j) = Ml(j) - U(i, j) * Ml(i);
        U(i, j) = U(i, j) + beta * Ml(j);
    end
    a = a - alpha * beta;
end
Lm = U.';
Um = U;
Dm = D;
end