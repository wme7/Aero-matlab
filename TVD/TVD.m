% Total Variation Diminishing (TVD) 
% by Manuel Diaz
% Example Algorithm to solve a 2D scalar advection equation:
%
% $$du/dt + dF(u)/dx + dG(u)/dy = 0$$
% where $dF/du = a$ and $dG/du = b$
%
clear; clc; close all;
%% Parameters
a = 0.5;  % Scalar velocity in x direction
b = 0.5;  % Scalar velocity in y direction
CFL  = 0.55;
sEnd = 10; % iterations

%% Domain Discretization
% Do nx and ny always > 4, for any 2D computation.
% Do nx or  ny == 1 to perform a 1D case.
dx = 1;     dy = 1; 
nx = 10;    ny = 1;
[x,y,d] = TVD_grid(nx,ny);

% Parameters of Blocks dimensions
l_1 = 1:ceil(nx/2); l_2 = ceil(nx/2)+1:nx;
h_1 = 1:ceil(ny/2); h_2 = ceil(ny/2)+1:ny;

%% Initial Condition: u(x,y,0)
u = zeros(ny,nx);
u(h_1,l_1) = 0.70; % block 1
%u(h_1,l_2) = 0.10; % block 2
%u(h_2,l_1) = 0.90; % block 3
u(h_2,l_2) = 0.50; % block 4

%% Constants
dt = min(dy,dx)*CFL/max(a,b);  % in this case dx = dy
nu_x = a*dt/dx;
nu_y = b*dt/dy;

%% Movie Parameters
%sEnd  = 500; % iterations
sPlot = 2; frames = sEnd/sPlot; counter = 0;
figure(1)
colordef white %black

%% TVD Method
time = 0;
for s = 1:sEnd
    % Compute Theta (smoothness coeficient)
    [theta_x,theta_y] = theta(u,nu_x,nu_y,d);

    % Flux Limiter
    [phi_x,phi_y] = van_leer(u,theta_x,theta_y,d);

    % Compute TVD Fluxes:
    [F_left,F_right,G_left,G_right]=TVD_flux(u,a,b,dx,dy,dt,phi_x,phi_y,d);
    
    % Compute new Step:
    u_next = u - dt/dx*(F_right-F_left) - dt/dy*(G_right-G_left);
    u = u_next;
    
    % Update Boundary Conditions: All Neumann BC's
    switch d % Dimension
        case{1} %1D problem
            u(y,1) = u(y,2); u(y,nx) = u(y,nx-1);
        case{2} %2D problem
            u(y,1) = u(y,2); u(y,nx) = u(y,nx-1);
            u(1,:) = u(2,:); u(ny,:) = u(ny-1,:);
    end
    
    % Update time
    time = time + dt;
    
    % Visualization
    if mod(s,sPlot) == 1
        counter = counter + 1;
        if d == 1
            plot(u)
        else
            contourf(u)
            colormap Autumn
            colorbar('location','southoutside')
        end
        M(counter) = getframe;
    end
    
    % Animated gif file
    if mod(s,sPlot) == 1
        F = getframe;
        if counter == 1
            [im,map] = rgb2ind(F.cdata,256,'nodither');
            im(1,1,1,sEnd/sPlot) = 0;
        end
        im(:,:,1,counter) = rgb2ind(F.cdata,map,'nodither');
    end
    
end

%% Simple Plot
figure(2)
if d == 1
    plot(u)
else
    contourf(u)
    colormap Autumn
    colorbar('location','southoutside')
end
title(['TVD Method, dx = ',num2str(dx),', dy = ',num2str(dy),', time: ',num2str(time)])
xlabel('x points')
ylabel('y points')

%% Make Movie
movie(M,2,10); % movie(M,n,fps)

%% Export to Gif
imwrite(im,map,'TVD.gif','DelayTime',0,'LoopCount',3)