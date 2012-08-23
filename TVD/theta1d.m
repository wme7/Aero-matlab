function [r_x] = theta1d(u,a)
%% Theta: smooth measuring factor (r_x)
% INPUTS
% u: Vector domain, u(j)
% a: x-velocity in our domain, a(j)

%% Compute r_x 
% Initialize variables
nx= length(u);
r_x = zeros(1,nx);

% Main loop
x = 2:nx-1;
for j = x
    % smooth measurement factor 'r_x'
    if u(j) == u(j+1)
        r_x(j) = 1;
    elseif a > 0 
        r_x(j) = (u(j) - u(j-1)) / (u(j+1) - u(j));
    elseif a < 0
        r_x(j) = (u(j+2) - u(j+1)) / (u(j+1) - u(j));
    end
end
r_x(1) = 1; r_x(nx) = 1;

return
        