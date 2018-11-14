function [Lm, Dm, Um] = modifyLDUlr(D, U, Ml, a)
[Lm, Dm, Um] = calcLDU(U.'*D*U + Ml*a*Ml');
end