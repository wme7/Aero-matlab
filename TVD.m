% Total Variation Diminishing (TVD) 
% Example Algorithm to solve a 2D scalar advection equation:

% du/dt + a*du/dx + b*du/dy = 0

clear; clc; close all;
%% Parameters
a = 0.5;  % Scalar velocity in x direction
b = 0.5;  % Scalar velocity in y direction
CFL  = 0.55;
%tEnd = 2.0; % out time
sEnd = 10; %iterations

%% Domain Discretization
dx = 0.1;     dy = 0.1;
nx = 50;    ny = 50;
x = 2:nx-1; y = 2:ny-1;

% Parameters of Blocks dimensions
d1 = 1:nx/2; d2 = nx/2+1:nx;
d3 = 1:ny/2; d4 = ny/2+1:ny;

%% Initial Condition: u(x,y,0)
u = zeros(ny,nx);
u(d3,d1) = 0.5;
u(d4,d2) = 0.7;

%% Constants
dt = min(dy,dx)*CFL/max(a,b);  % in this case dx = dy
nu_x = a*dt/dx;
nu_y = b*dt/dy;
vxp = max(a,0);
vxm = min(a,0);
vyp = max(b,0);
vym = min(b,0);

%% Movie Parameters
%sEnd  = 500; % iterations
sPlot = 2; frames = sEnd/sPlot; counter = 0;
figure(1)
colordef white %black

%% TVD Method
time = 0;
%while time < tEnd
for s = 1:sEnd
    % Compute Theta
    theta_x = zeros(ny,nx);
    theta_y = zeros(ny,nx);
    for  j = y
        for i = x
            if u(j,i+1) == u(j,i)
                theta_x(j,i) = 0;
            else
                theta_x(j,i) = (u(j,i+1-sign(nu_x)) - u(j,i-sign(nu_x))) ...
                    /(u(j,i+1) - u(j,i));
            end
            if u(j+1,i) == u(j,i)
                theta_y(j,i) = 0;
            else
                theta_y(j,i) = (u(j+1-sign(nu_y),i) - u(j-sign(nu_y),i)) ...
                    /(u(j+1,i) - u(j,i));
            end
        end
    end

    % Flux Limiter : Van Leer
    phi_x = zeros(ny,nx);
    phi_y = zeros(ny,nx);
    for j = y
        for i = x
            if theta_x == 0 
                phi_x(j,i) = 0;
            else
                phi_x(j,i) = (abs(theta_x(j,i)) + theta_x(j,i)) ... 
                    / (1 + abs(theta_y(j,i)));
            end 
            if theta_y == 0
                phi_y(j,i) = 0;
            else
                phi_y(j,i) = (abs(theta_y(j,i)) + theta_y(j,i)) ...
                    / (1 + abs(theta_y(j,i)));
            end 
        end
    end

    %Compute TVD Fluxes:
    fsl = zeros(ny,nx); fsh = zeros(ny,nx); fs = zeros(nx,ny);
    frl = zeros(ny,nx); frh = zeros(ny,nx); fr = zeros(nx,ny);
    gsl = zeros(ny,nx); gsh = zeros(ny,nx); gs = zeros(nx,ny);
    grl = zeros(ny,nx); grh = zeros(ny,nx); gr = zeros(nx,ny);
    for j = y
        for i = x
            % FLux to the left (s)
            fsl(i,j) = vxp*u(i-1,j) + vxm*u(i,j);
            fsh(i,j) = 1/2*a*(u(i-1,j)+u(i,j)) ...
                - 1/2*dt/dx*(a^2)*(u(i,j)-u(i-1,j)); 
            fs(i,j)  = fsl(i,j) + (phi_x(i-1,j)*(fsh(i,j) - fsl(i,j)));
            % Flux to the right (r)
            frl(i,j) = vxp*u(i,j) + vxm*u(i+1,j);
            frh(i,j) = 1/2*a*(u(i,j)+u(i+1,j)) ...
                - 1/2*dt/dx*(a^2)*(u(i+1,j)-u(i,j));
            fr(i,j)  = frl(i,j) + (phi_x(i,j)*(frh(i,j) - frl(i,j)));
            % FLux to the up (s)
            gsl(i,j) = vyp*u(i,j-1) + vym*u(i,j);
            gsh(i,j) = 1/2*b*(u(i,j-1)+u(i,j)) ...
                -1/2*dt/dy*(b^2)*(u(i,j)-u(i,j-1));
            gs(i,j)  = gsl(i,j) + (phi_y(i,j-1)*(gsh(i,j) - gsl(i,j)));
            % FLux to the down (s)
            grl(i,j) = vyp*u(i,j) + vym*u(i,j+1);
            grh(i,j) = 1/2*b*(u(i,j)+u(i,j+1)) ...
                -1/2*dt/dy*(b^2)*(u(i,j+1)-u(i,j));
            gr(i,j)  = grl(i,j) + (phi_y(i,j)*(grh(i,j) - grl(i,j)));
        end
    end 
    u_next = u - dt/dx*(fr-fs) - dt/dy*(gr-gs);
    u = u_next;
    % Update Boundary Conditions: All Neumann BCs
    u(y,1) = u(y,2); u(y,nx) = u(y,nx-1);
    u(1,:) = u(2,:); u(ny,:) = u(ny-1,:);
    
    % Update time
    time = time + dt;
    
    % Visualization
    if mod(s,sPlot) == 1
        counter = counter + 1;
        contourf(u)
        colormap Autumn
        colorbar('location','southoutside')
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
contourf(u)
title(['TVD Method, dx = ',num2str(dx),', dy = ',num2str(dy),', time: ',num2str(time)])
xlabel('x points')
ylabel('y points')
colormap Autumn
colorbar('location','southoutside')

%% Make Movie
%movie(M,2,10); % movie(M,n,fps)

%% Export to Gif
%imwrite(im,map,'TVD.gif','DelayTime',0,'LoopCount',3)