function [x, niters] = PCG(A, b, x0)

% Perform incomplete Cholesky factorization on Asparse
L = ichol(sparse(A), struct('type','ict','droptol',1e-4,'michol','off'));
M = L*L';

% Conjugate Gradient Method without preconditioning
x = x0;
r = b; 
niters = 0;
zold = M\r;


while norm(r) > eps()*norm(b) && niters < 1000
    z = M\r;
    if niters == 0
        p = z;
    else    
        gamma = dot(r',z)/dot(rold',zold);
        p = z + gamma*p;
    end
    q = A*p;
    alpha = dot(r',z)/dot(p',q);
    x = x + alpha * p; 
    rold = r;
    zold = z;
    r = r - alpha * q; 
    niters = niters + 1; 
end
end