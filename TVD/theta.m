function [theta_x,theta_y] = theta(u,nu_x,nu_y,Dimension)
% Compute Theta
[ny,nx] = size(u);
theta_x = zeros(ny,nx);
theta_y = zeros(ny,nx);
switch Dimension
    
    case{1} % 1D Problem
        x = 2:nx-1;
        for i = x
           if u(i+1) == u(i)
               theta_x(i) = 0;
           else
               theta_x(i) = (u(i+1-sign(nu_x)) - u(i-sign(nu_x))) ...
                   /(u(i+1) - u(i));
           end
        end
        %theta_d = theta_x;
        
    case{2} % 2D Problem
        x = 2:nx-1;
        y = 2:ny-1;
        for  j = y
            for i = x
                if u(j,i+1) == u(j,i)
                    theta_x(j,i) = 0;
                else
                    theta_x(j,i) = (u(j,i+1-sign(nu_x)) - u(j,i-sign(nu_x))) ...
                        /(u(j,i+1) - u(j,i));
                end
                if u(j+1,i) == u(j,i)
                    theta_y(j,i) = 0;
                else
                    theta_y(j,i) = (u(j+1-sign(nu_y),i) - u(j-sign(nu_y),i)) ...
                        /(u(j+1,i) - u(j,i));
                end
            end
        end
        %theta_d = [theta_x theta_y];
    otherwise
            error('only supported cases d=1 and d=2');
end
return
        