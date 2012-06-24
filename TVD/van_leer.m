function [phi_x,phi_y] = van_leer(u,theta_x,theta_y,Dimension)
[ny,nx] = size(u);
phi_x = zeros(ny,nx);
phi_y = zeros(ny,nx);
switch Dimension
    
    case{1} % 1D Problem
        x = 2:nx-1;
        for i = x
            if theta_x == 0 
                phi_x(i) = 0;
            else
                phi_x(i) = (abs(theta_x(i)) + theta_x(i)) ... 
                    / (1 + abs(theta_y(i)));
            end
        end
                
    case{2} % 2D Problem
        x = 2:nx-1;
        y = 2:ny-1;
        for j = y
            for i = x
                if theta_x == 0 
                    phi_x(j,i) = 0;
                else
                    phi_x(j,i) = (abs(theta_x(j,i)) + theta_x(j,i)) ... 
                        / (1 + abs(theta_y(j,i)));
                end 
                if theta_y == 0
                    phi_y(j,i) = 0;
                else
                    phi_y(j,i) = (abs(theta_y(j,i)) + theta_y(j,i)) ...
                        / (1 + abs(theta_y(j,i)));
                end 
            end
        end
        
    otherwise
            error('only supported cases d=1 and d=2');
end
return