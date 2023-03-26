function [x, niters] = CG(A, b, x0)
% Conjugate Gradient Method without preconditioning
x = x0;
r = b; 
p = r; 
niters = 0;
rold = r;

while norm(r) > eps()*norm(b) && niters < 1000
    if niters ~= 0
        gamma = dot(r',r)/dot(rold',rold);
        p = r + gamma*p;
    end    
    alpha = (r' * r) / (p' * A * p);
    x = x + alpha * p; 
    rold = r; 
    r = r - alpha * A * p; 
    niters = niters + 1; 
end
end