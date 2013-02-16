function [r_x] = theta1d(u,a)
%% Theta: smooth measuring factor (r_x)
% INPUTS
% u: Vector domain, u(j)
% a: x-velocity in our domain, a(j)

%% Compute r_x 
% Initialize variables
nx  = length(u);
r_x = zeros(1,nx);

%% Test a 
% Check whether a is scalar or a vector of velocities with
% prescribed velocities in the entire domain
[m,n,r] = size(a);
if m == 1 && n == 1 && r == 1
    a = a*ones(1,nx); % map the x-velocity 
else
    %do nothing
end

% Main loop
x = 2:nx-1;
for j = x
    % smooth measurement factor 'r_x'
    if u(j) == u(j+1)
        r_x(j) = 1;
    elseif a(j) > 0 
        r_x(j) = (u(j) - u(j-1)) / (u(j+1) - u(j));
    elseif a(j) < 0
        u(nx+1) = u(nx); % we will need an extra column value
        r_x(j) = (u(j+2) - u(j+1)) / (u(j+1) - u(j));
    end
end
r_x(1) = 1; r_x(nx) = 1;

return