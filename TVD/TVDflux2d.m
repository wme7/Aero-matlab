function [f_s,f_r,g_s,g_r] = TVDflux2d(u,a,b,dtdx,dtdy,phi_x,phi_y)
%% Grid Size
[ny,nx] = size(u);
x = 2:nx-1;
y = 2:ny-1;

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

% Compute Flux Spliting Factors
a_p = max(a,0); % a{+}
a_m = min(a,0); % a{-}
b_p = max(b,0); % b{+}
b_m = min(b,0); % b{-}

%% Compute TVD fluxes
% Note:
% s = Left flux & r = Right flux
% l = low flux  & h = high flux

% Initalize Arrays
fsl = zeros(ny,nx); fsh = zeros(ny,nx); f_s = zeros(ny,nx);
frl = zeros(ny,nx); frh = zeros(ny,nx); f_r = zeros(ny,nx);
gsl = zeros(ny,nx); gsh = zeros(ny,nx); g_s = zeros(ny,nx);
grl = zeros(ny,nx); grh = zeros(ny,nx); g_r = zeros(ny,nx);

    
for j = y
    for i = x
        % Flux to the right(r)
        frl(j,i) = a_p(j,i)*u(j,i) + a_m(j,i)*u(j,i+1);
        frh(j,i) = 1/2*a(j,i)*(u(j,i) + u(j,i+1)) ...
            - 1/2*dtdx*(a(j,i)^2)*(u(j,i+1)-u(j,i));
        f_r(j,i) = frl(j,i) + phi_x(j,i)*(frh(j,i) - frl(j,i));
        
        % FLux to the left (s)
        fsl(j,i) = a_p(j,i)*u(j,i-1) + a_m(j,i)*u(j,i);
        fsh(j,i) = 1/2*a(j,i)*(u(j,i-1) + u(j,i)) ...
            - 1/2*dtdx*(a(j,i)^2)*(u(j,i)-u(j,i-1)); 
        f_s(j,i) = fsl(j,i) + phi_x(j,i-1)*(fsh(j,i) - fsl(j,i));
            
        % FLux to the up (r)
        grl(j,i) = b_p(j,i)*u(j,i) + b_m(j,i)*u(j+1,i);
        grh(j,i) = 1/2*b(j,i)*(u(j,i) + u(j+1,i)) ...
            - 1/2*dtdy*(b(j,i)^2)*(u(j+1,i)-u(j,i));
        g_r(j,i) = grl(j,i) + phi_y(j,i)*(grh(j,i) - grl(j,i));

        % FLux to the down (s)
        gsl(j,i) = b_p(j,i)*u(j-1,i) + b_m(j,i)*u(j,i);
        gsh(j,i) = 1/2*b(j,i)*(u(j-1,i) + u(j,i)) ...
            - 1/2*dtdy*(b(j,i)^2)*(u(j,i)-u(j-1,i));
        g_s(j,i) = gsl(j,i) + phi_y(j-1,i)*(gsh(j,i) - gsl(j,i));
        
    end
end

return