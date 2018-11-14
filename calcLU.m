function [L, U] = calcLU(A)
N = size(A,1); 
for p=1:N-1
    A(p, p+1:end) = A(p, p+1:end) / A(p, p);
    for i=find(A(p+1:end, p))+p
        A(i, p+1:end) = A(i, p+1:end) - A(i, p) * A(p, p+1:end);
    end
end
L = tril(A);
U = triu(A, 1) + speye(N);
