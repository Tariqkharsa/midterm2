function [x, niters] = Method_of_Steepest_Descent(A, b, x0)
% Solve the linear system Ax=b using the Method of Steepest Descent.

% Initialize variables
x = x0;
r = b - A*x;
niters = 0;

% Perform iterations until convergence or maximum number of iterations
while norm(r) > eps()*norm(b) && niters < 1000
    niters = niters + 1;
    q = A*r;
    alpha = dot(r', r)/dot(r', q);
    x = x + alpha*r;
    r = r - alpha*q;
end

end


