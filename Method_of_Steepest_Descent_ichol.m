function [x, niters] = Method_of_Steepest_Descent_ichol(A, b, x0)

% Perform incomplete Cholesky factorization on Asparse
L = ichol(sparse(A), struct('type','ict','droptol',1e-3,'michol','off'));
M = L*L';

% Set up initial guess and residual
x = x0;
r = b - A*x;
niters = 0;

while norm(r) > eps()*norm(b) && niters < 1000
    p = M\r;
    q = A*p;
    alpha = dot(p', r) / dot(p', q);
    x = x + alpha*p;
    r = r - alpha*q;
    niters = niters + 1;
end

end