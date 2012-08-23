%% 2D Total Variation Diminishing (TVD) subroutine
% by Manuel Diaz
%
% to solve a 2D scalar advection equation:
%
% $$du/dt + dF(u)/dx + dG(u)/dy = 0$$
% where $dF/du = a$ and $dG/du = b$
%
clear; clc; close all;

%% Main Parameters
    a = 0.45;  % Scalar velocity in x direction
    b = 0.65;  % Scalar velocity in y direction
  cfl = 0.60;  % CFL condition
t_end = 0.50;  % iterations

%% Domain 
    d = 2; % 2D domain is used
% 2D domain where dx = dy, nx > 1 and ny > 4!!
   nx = 10;   ny = 4;
[x,dx,y,dy] = grid2d(0,9,nx,0,3,ny);

% time discreti zation
   dt = min(dy,dx)*cfl/max(a,b); 
    t = 0:dt:t_end;
 dtdx = dt/dx;
 dtdy = dt/dy;

%% Initial Condition: u(x,y,0)
% Domain Velocity
a = a*ones(ny,nx); % map the x-velocity 
b = b*ones(ny,nx); % map the y-velocity 
a_p = max(0,a); a_m = min(0,a);
b_p = max(0,b); b_m = min(0,b);

% Parameters of regions dimensions
x_middle = ceil(nx/2);
y_middle = ceil(ny/2);
l_1 = 1:x_middle; l_2 = x_middle+1:nx;
h_1 = 1:y_middle; h_2 = y_middle+1:ny;

% Initial Condition for our 2D domain
u = zeros(ny,nx);
u_0 = zeros(ny,nx);
u_0(h_1,l_1) = 0.70; % region 1
u_0(h_1,l_2) = 0.10; % region 2
u_0(h_2,l_1) = 0.90; % region 3
u_0(h_2,l_2) = 0.50; % region 4

% %% Movie Parameters
% %sEnd  = 500; % iterations
% sPlot = 2; frames = sEnd/sPlot; counter = 0;
% figure(1)
% colordef white %black

%% Main Loop 
% load Initial condition into our Domain
u = u_0; 

% Initialize Matrix Arrays
u_next(ny,nx);

r_x = zeros(ny,nx);
r_y = zeros(ny,nx);

F_rl = zeros(ny,nx);
F_rh = zeros(ny,nx);
F_ll = zeros(ny,nx);
F_lh = zeros(ny,nx);
F_right = zeros(ny,nx);
F_left  = zeros(ny,nx);

G_rl = zeros(ny,nx);
G_rh = zeros(ny,nx);
G_ll = zeros(ny,nx);
G_lh = zeros(ny,nx);
G_right = zeros(ny,nx);
G_left  = zeros(ny,nx);

for k = t
    % Compute Theta (smoothness coeficient)
    [r_x,r_y] = theta(u,nu_x,nu_y,d);

    % Compute flux Limiter 
    [phi_x,phi_y] = van_leer(u,r_x,r_y,d);

    % Compute TVD Fluxes:
    [F_left,F_right,G_left,G_right]=TVD_flux(u,a,b,dx,dy,dt,phi_x,phi_y,d);
    
    % Compute new Step:
    u_next = u - dt/dx*(F_right-F_left) - dt/dy*(G_right-G_left);
    
    % Update Information
    u = u_next;
    
    % Update Boundary Conditions: All Neumann BC's
   	u(y,1) = u(y,2); u(y,nx) = u(y,nx-1);
    u(1,:) = u(2,:); u(ny,:) = u(ny-1,:);
    
    % Update Iteration counter
    % s = s + 1;
    
%     % Visualization
%     if mod(s,sPlot) == 1
%         counter = counter + 1;
%         if d == 1
%             plot(u)
%         else
%             contourf(u)
%             colormap Autumn
%             colorbar('location','southoutside')
%         end
%         M(counter) = getframe;
%     end
%     
%     % Animated gif file
%     if mod(s,sPlot) == 1
%         F = getframe;
%         if counter == 1
%             [im,map] = rgb2ind(F.cdata,256,'nodither');
%             im(1,1,1,sEnd/sPlot) = 0;
%         end
%         im(:,:,1,counter) = rgb2ind(F.cdata,map,'nodither');
%     end
    
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

% %% Make Movie
% movie(M,2,10); % movie(M,n,fps)
% 
% %% Export to Gif
% imwrite(im,map,'TVD.gif','DelayTime',0,'LoopCount',3)