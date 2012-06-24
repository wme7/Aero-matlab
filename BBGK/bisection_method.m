function x = bisection_method(f,a,b,tol,Nmax)
%% bisection method 
% Numerical Implementation of a bisection method:
g = f(x);
N = 1
while N <= Nmax
    c = (a+b)/2; % Midpoint
    if g(c) = 0 or (b-a)/2 < tol % Solution found
        x = c;
        break
    else
        N = N + 1;
        
    
end

