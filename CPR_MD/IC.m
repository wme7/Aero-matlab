function u0 = IC(x,ICcase)
% Function to construct the Initial Condition

% Create the selected IC
switch ICcase
    case {1} % Gaussian wave
        u0 = exp(-(6*(x-0.5)).^2);
        
    case {2} % sinusoidal wave
        %u0 = 1/2 + sin(2*pi*x);
        u0 = 2 + 0.5*sin(2*pi*x);
    case {3} % Riemann problem
             % u = 1 for x <  x_mid
             % u = 0 for x >= x_mid
        u0 = zeros(size(x));
        rhs = find(x<x_mid);
        u0(rhs) = 1;
    otherwise
        error('case not in the list')
end