function Sm = generateSm(idxM, U)
N = size(U, 1);
S = zeros(N, 1);
n = 0;
for i = 1:size(idxM, 1)
    p = idxM(i);
    t = [p p];
    while (size(t, 2) > 1 && ~ismember(t(2), S))
        n = n + 1;
        S(n) = t(2);
        p = t(2);
        t = find(U(p, :), 2);
    end
end
Sm = S(1:n);
end