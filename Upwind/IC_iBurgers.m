function u0 = IC_iBugers(x,ICcase)
% Function to construct the Initial Condition for any FDM

%% Calculate the number of points used in the computation.
npt = length(x);

%% Parameters
x_min = x(1);   % x value at the begining of our domain
x_max = x(end); % x value at the end of our domain

x_diff  = (x_max-x_min);        % length of spatial range
x_mid   = (x_min+x_max)/2.0;    % mid-point of spatial range

v_min = +0.20;  % maximum velocity
v_max = +0.70;  % minimum velocity

f_tip   = 0.30; % location of gaussian and triangle
f_gauss = 0.40; % width parameter for gaussian
f_width = 0.20; % width parameter for triangle and slope

%% Create the selected
switch ICcase
    
    case {1} % Gaussian
        for i=1:npt
            x_tip   = x_min + f_tip*x_diff;  % center of triangle
            u0(1,i) = (v_max-v_min) * exp( -( (x(i)-x_tip) / f_gauss*x_diff )^2 ) + v_min;
            
        end
                
    case {2} % Slope
        for i=1:npt
            x_left  = x_mid-f_width*x_diff;
            x_right = x_mid+f_width*x_diff;
            u0(1,i) = v_max;
            if x(i) > (x_left)
                if x(i) < (x_right)
                    u0(1,i) = v_max - (v_max-v_min)*((x(i)-x_left)/(x_right-x_left));
                else
                    u0(1,i) = v_min;
                end
            end
        end
        
    case {3} % Triangle
        for i=1:npt
            x_tip   = x_min + f_tip*x_diff;  % center of triangle
            x_left  = x_tip-f_width*x_diff;  % left border of triangle
            x_right = x_tip+f_width*x_diff;  % right border of triangle
            u0(1,i) = v_min;                  % velocity outside triangle region
            if x(i) > x_left
                if x(i) < x_right
                    if x(i) < x_tip              % left side of triangle
                        u0(1,i) = v_min + (v_max-v_min)*((x(i)-x_left)/(x_tip-x_left));
                    else                         % right side of triangle
                        u0(1,i) = v_max - (v_max-v_min)*((x(i)-x_tip)/(x_right-x_tip));
                    end
                end
            end
        end

    case {4} % sinusoidal wave
        
        %u0 = 1/2 + sin(2*pi*x);
        u0 = 2 + 0.5*sin(pi*x);
        
    otherwise
        error('case not in the list')
end