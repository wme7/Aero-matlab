function u0 = IC(x,ICcase)
% Function to construct the Initial Condition

% Create the selected IC
switch ICcase
    case 1 % Gaussian wave
        xmid = (x(end)-x(1))/2;
        u0 = exp(-(2*(x-xmid)).^2);
    case 2 % sinusoidal wave
        u0 = 0.5 + sin(x);
    case 3 % Riemann problem
             % u = 1 for x <  x_mid
             % u = 0 for x >= x_mid
        u0 = ones(size(x));
        xmid = 0.45;%(x(end)-x(1))/2;
        rhs = find(x<=xmid);
        u0(rhs) = 2; %1.2;
    case 4 % centered sinusoidal wave
        u0 = sin(x);
    case 5 % Tanh
        a = x(1); b = x(end);
        xi = (4-(-4))/(b-a)*(x - a) - 4;
        u0 = 1/2*(tanh(-4*xi)+1);
    otherwise
        error('case not in the list')
end