function [r_x,r_y] = theta2d(u,a,b)
%% Theta: smooth measuring factor (r)
% INPUTS
% u: matrix domain, u(ny,nx)
% a: x-velocity in domain, a(ny,nx)
% b: y-velocity in domain, b(ny,nx)

%% Grid zise
[ny,nx] = size(u);
r_x = zeros(ny,nx);
r_y = zeros(ny,nx);

%% Test if a and b 
% Check whether a and b are scalar velocities or 
% prescribed velocities in the entire domain
[m,n,r] = size(a);
if m == 1 && n == 1 && r == 1
    a = a*ones(ny,nx); % map the x-velocity 
else
    %do nothing
end
clear m n r;

[m,n,r] = size(b);
if m == 1 && n == 1 && r == 1
    b = b*ones(ny,nx); % map the y-velocity 
else
    %do nothing
end
clear m n r;


%% Compute r_x and r_y
x = 2:nx-1;
for j = 1:ny
    for i = x
        % smooth measurement factor 'r_x'
        if u(j,i) == u(j,i+1)
            r_x(j,i) = 1;
        elseif a(j,i) > 0
            r_x(j,i) = (u(j,i) - u(j,i-1))/(u(j,i+1) - u(j,i));
        elseif a(j,i) < 0
            u(j,nx+1) = u(j,nx); % we will need an extra column value
            r_x(j,i) = (u(j,i+2) - u(j,i+1))/(u(j,i+1) - u(j,i));
        end
    end
end
r_x(:,1) = 1; r_x(:,nx) = 1;

y = 2:ny-1;
for j = y
    for i = 1:nx
        % smooth measurement factor 'r_y'
        if u(j,i) == u(j+1,i)
            r_y(j,i) = 1;
        elseif b(j,i) > 0
            r_y(j,i) = (u(j,i) - u(j-1,i))/(u(j+1,i) - u(j,i));
        elseif b(j,i) < 0
            u(ny+1,i) = u(ny,i); % we will need an extra row value
            r_y(j,i) = (u(j+2,i) - u(j+1,i))/(u(j+1,i) - u(j,i));
        end
    end
end
r_y(1,:) = 1; r_y(ny,:) = 1;

return     