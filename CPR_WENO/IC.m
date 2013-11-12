function u0 = IC(x,ICcase)
% Function to construct the Initial Condition

% Create the selected IC
switch ICcase
    case {1} % Gaussian wave
        u0 = exp(-(6*(x-0.5)).^2);
        
    case {2} % sinusoidal wave
        u0 = 0.5*(2 + sin(2*pi*x));
    case {3} % Riemann problem
             % u = 1 for x <  x_mid
             % u = 0 for x >= x_mid
        u0 = ones(size(x));
        xmid = (x(end)-x(1))/2;
        rhs = find(x<=xmid);
        u0(rhs) = 2;
    otherwise
        error('case not in the list')
end