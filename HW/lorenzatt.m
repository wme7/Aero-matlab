% Lorenz Attractor equations solved by ODE Solve
%% x' = sigma*(y-x)
%% y' = x*(rho - z) - y
%% z' = x*y - beta*z
function dx = lorenzatt(X)
    rho = 28; sigma = 10; beta = 8/3;
    dx = zeros(3,1);
    dx(1) = sigma*(X(2) - X(1));
    dx(2) = X(1)*(rho - X(3)) - X(2);
    dx(3) = X(1)*X(2) - beta*X(3);
    return
end