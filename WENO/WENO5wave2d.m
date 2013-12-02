%% Weighted Essentially non-Oscilaroty (WENO) Routine
% to solve the scalar advection equation:
%
% $$du/dt + dF(u)/dx + dG(u)/dy = 0$$
%
% where $dF/du = a$ and $dG/du = b$
%
% by Manuel Diaz, manuel.ade'at'gmail.com 
% Institute of Applied Mechanics, 2012.08.21
%
clear; clc; close all;
%% Main Parameters
      a = 0.40; % Scalar velocity in x direction
      b =-0.60; % Scalar velocity in y direction
    cfl = 0.40; % CFL condition
  t_end = 2.50; % Final time
      k = 3;    % WENO Order: 1st, 2nd and 3rd Orders available.  

%% Domain 
    d = 2; % 2D domain is used
% 2D domain where dx = dy, nx > 1 and ny > 4!!
   nx = 40;   ny = 40;
[x,dx,y,dy] = grid2d(0,7,nx,0,15,ny);

% time discretization
   dt = min(dy,dx)*cfl/max(abs(a),abs(b)); 
    t = 0:dt:t_end; 
 dtdx = dt/dx; % precomputed to save flops
 dtdy = dt/dy; % precomputed to save flops

% Domain Velocity
a_p = max(0,a); a_m = min(0,a);
b_p = max(0,b); b_m = min(0,b);

%% Initial Condition (IC)
% Build IC
u_0 = u_zero2d(x,y);

%% Movie Parameters
sEnd  = length(t); % iterations
sPlot = 1; frames = sEnd/sPlot; counter = 0;
figure(1)
colordef white %black

%% Main Loop 
% iteration counter
s = 1;

% load Initial condition into our Domain
u = u_0; 

% Initialize Matrix Arrays
u_next = zeros(ny,nx);

for kk = t
    % Compute WENO Fluxes:
    [F_l,F_r,G_l,G_r] = WENOflux2d(u,a,b); % is underconstruction!
    
    % Compute new Step:
    u_next = u - dtdx*(F_r - F_l) - dtdy*(G_r - G_l);
    
    % Update BCs: All Neumann BC's
   	u_next(:,1)  = u_next(:,4); 
    u_next(:,2)  = u_next(:,4);
    u_next(:,3)  = u_next(:,4);
    
    u_next(:,nx  ) = u_next(:,nx-3);
    u_next(:,nx-1) = u_next(:,nx-3);
    u_next(:,nx-2) = u_next(:,nx-3);
    
    u_next(1,:)  = u_next(4,:);
    u_next(2,:)  = u_next(4,:);
    u_next(3,:)  = u_next(4,:);
    
    u_next(ny  ,:) = u_next(ny-3,:);
    u_next(ny-1,:) = u_next(ny-3,:);
    u_next(ny-2,:) = u_next(ny-3,:);
    
    % Update Information
    u = u_next;
    
    % Update Iteration counter
    s = s + 1;
    
    % Visualization
    if mod(s,sPlot) == 0 % I want 1 frame per iteration
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
    if mod(s,sPlot) == 0 % I want 1 frame per iteration
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


%% Make Movie
movie(M,2,10); % movie(M,n,fps)

% %% Export to Gif
imwrite(im,map,'WENO3_2D.gif','DelayTime',0,'LoopCount',5)