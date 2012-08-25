%% 2D Total Variation Diminishing (TVD) 2D subroutine
% to solve the scalar advection equation:
%
% $$du/dt + dF(u)/dx + dG(u)/dy = 0$$
%
% where $dF/du = a$ and $dG/du = b$
%
% Using flux limiter functions: 
% Cases: {1} Van Leer
%        {2} Superbee
%        {3} Minmod
%        {4} Koren
%
% by Manuel Diaz, manuel.ade'at'gmail.com 
% Institute of Applied Mechanics, 2012.08.25

clear; clc; close all;

%% Main Parameters
      a = 0.60; % Scalar velocity in x direction
      b =-0.40; % Scalar velocity in y direction
    cfl = 0.60; % CFL condition
  t_end = 2.00; % iterations
limiter = 2;    % Options: 1(Vl), 2(Sb), 3(Mm), 4(koren)

%% Domain 
    d = 2; % 2D domain is used
% 2D domain where dx = dy, nx > 1 and ny > 4!!
   nx = 8;   ny = 16;
[x,dx,y,dy] = grid2d(0,7,nx,0,15,ny);

% time discretization
   dt = min(dy,dx)*cfl/max(abs(a),abs(b)); 
    t = 0:dt:t_end; 
 dtdx = dt/dx; % precomputed to save flops
 dtdy = dt/dy; % precomputed to save flops

%% Initial Condition: u(x,y,0)
% Domain Velocity
a_p = max(0,a); a_m = min(0,a);
b_p = max(0,b); b_m = min(0,b);

% Parameters of regions dimensions
x_middle = ceil(nx/2);
y_middle = ceil(ny/2);
l_1 = 1:x_middle; l_2 = x_middle+1:nx;
h_1 = 1:y_middle; h_2 = y_middle+1:ny;

% Initial Condition for our 2D domain
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
u_next = zeros(ny,nx);

rx = zeros(ny,nx);
ry = zeros(ny,nx);

F_r = zeros(ny,nx);
F_l = zeros(ny,nx);

G_r = zeros(ny,nx);
G_l = zeros(ny,nx);

for k = t
    % Compute Theta (smoothness coeficient)
    [rx,ry] = theta2d(u,a,b);
    
    % Compute flux Limiter 
    [phix,phiy] = fluxlimiter2d(rx,ry,limiter);
        
    % Compute TVD Fluxes:
    [F_l,F_r,G_l,G_r] = flux2d(u,a,b,dtdx,dtdy,phix,phiy);
    
    % Compute new Step:
    u_next = u - dtdx*(F_r - F_l) - dtdy*(G_r - G_l);
    
    % Update BCs: All Neumann BC's
   	u_next(:,1)  = u_next(:,2); 
    u_next(:,nx) = u_next(:,nx-1);
    u_next(1,:)  = u_next(2,:); 
    u_next(ny,:) = u_next(ny-1,:);
    
    % Update Information
    u = u_next;
    
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
subplot(1,2,1)
hold on
    %contourf(u)
    surface(u)
    colormap Autumn
    colorbar('location','southoutside')
title(['TVD Method, dx = ',num2str(dx),', dy = ',num2str(dy),', time: ',num2str(t_end)])
xlabel('x points'); ylabel('y points')
hold off

subplot(1,2,2)
hold on
    %contourf(u_0)
    surface(u_0)
    colormap Autumn
    colorbar('location','southoutside')
title(['Initial Condition, dx = ',num2str(dx),', dy = ',num2str(dy),', time: ',num2str(t(1))])
xlabel('x points'); ylabel('y points')
hold off


% %% Make Movie
% movie(M,2,10); % movie(M,n,fps)
% 
% %% Export to Gif
% imwrite(im,map,'TVD.gif','DelayTime',0,'LoopCount',3)