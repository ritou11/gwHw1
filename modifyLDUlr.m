function [Lm, Dm, Um] = modifyLDUlr(L, D, U, Ml, a)
[Lm, Dm, Um] = calcLDU(L*D*U + Ml*a*Ml');
end