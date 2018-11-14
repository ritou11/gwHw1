function [Lm, Dm, Um] = modifyLDUr1(L, D, U, Ml, a)
[Lm, Dm, Um] = calcLDU(L*D*U + Ml*a*Ml');
end