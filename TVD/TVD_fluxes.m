function [f_s,f_r,g_s,g_r] = TVD_flux(u,a,b,dx,dy,dt,phi_x,phi_y,Dimension)
%% Compute TVD fluxes

% Grid Size
[ny,nx] = size(u);
x = 2:nx-1;
y = 2:ny-1;

% Flux Spliting Factors
vxp = max(a,0);
vxm = min(a,0);
vyp = max(b,0);
vym = min(b,0);

% Initalize Arrays
fsl = zeros(ny,nx); fsh = zeros(ny,nx); f_s = zeros(ny,nx);
frl = zeros(ny,nx); frh = zeros(ny,nx); f_r = zeros(ny,nx);
gsl = zeros(ny,nx); gsh = zeros(ny,nx); g_s = zeros(ny,nx);
grl = zeros(ny,nx); grh = zeros(ny,nx); g_r = zeros(ny,nx);

switch Dimension
    
    case{1} % 1D Problem
        
        for i = x
            % FLux to the left (s)
            fsl(i) = vxp*u(i-1) + vxm*u(i);
            fsh(i) = 1/2*a*(u(i-1) + u(i)) ...
                - 1/2*dt/dx*(a^2)*(u(i)-u(i-1)); 
            f_s(i)  = fsl(i) + (phi_x(i-1)*(fsh(i) - fsl(i)));
           
            % Flux to the right (r)
            frl(i) = vxp*u(i) + vxm*u(i+1);
            frh(i) = 1/2*a*(u(i) + u(i+1)) ...
                - 1/2*dt/dx*(a^2)*(u(i+1)-u(i));
            f_r(i) = frl(i) + (phi_x(i)*(frh(i) - frl(i)));
        end
        
    case{2} % 2D Problem    
        
        for j = y
            for i = x
            % FLux to the left (s)
            fsl(j,i) = vxp*u(j,i-1) + vxm*u(j,i);
            fsh(j,i) = 1/2*a*(u(j,i-1) + u(j,i)) ...
                - 1/2*dt/dx*(a^2)*(u(j,i)-u(j,i-1)); 
            f_s(j,i)  = fsl(j,i) + (phi_x(j,i-1)*(fsh(j,i) - fsl(j,i)));
            
            % Flux to the right (r)
            frl(j,i) = vxp*u(j,i) + vxm*u(j,i+1);
            frh(j,i) = 1/2*a*(u(j,i) + u(j,i+1)) ...
                - 1/2*dt/dx*(a^2)*(u(j,i+1)-u(j,i));
            f_r(j,i) = frl(j,i) + (phi_x(j,i)*(frh(j,i) - frl(j,i)));
            
            % FLux to the down (s)
            gsl(j,i) = vyp*u(j-1,i) + vym*u(j,i);
            gsh(j,i) = 1/2*b*(u(j-1,i) + u(j,i)) ...
                -1/2*dt/dy*(b^2)*(u(j,i)-u(j-1,i));
            g_s(j,i) = gsl(j,i) + (phi_y(j-1,i)*(gsh(j,i) - gsl(j,i)));
            
            % FLux to the up (r)
            grl(j,i) = vyp*u(j,i) + vym*u(j+1,i);
            grh(j,i) = 1/2*b*(u(j,i) + u(j+1,i)) ...
                -1/2*dt/dy*(b^2)*(u(j+1,i)-u(j,i));
            g_r(j,i) = grl(j,i) + (phi_y(j,i)*(grh(j,i) - grl(j,i)));
            end
        end
    otherwise 
        error('Only supported case d = 1 and d = 2')
end
return
