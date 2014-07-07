%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Solving 2-D heat equation with Nodal DG
%
%      du/dt = k*( d^2u/dx^2 + d^2u/dy^2 ),  in (x,y) \in [a,b]x[c,d]
%                 and where u = u(x,y)
%
%              coded by Manuel Diaz, NTU, 2013.07.30
%                               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;

%% Preprocessing

% Parameters
 cfl = 0.2;
   N = 50;
tEnd = 0.05;

% build grid
x = linspace(-1,1,N+1); dx = 1/N;
y = linspace(-1,1,N+1); dy = 1/N;
[x,y] = meshgrid(x,y);

% set thermal diffusivity coef
k=0.1; 

% Set IC
u0 = IC2d(x,y,2);

% set BCs
u0(1,1:N) = 0.1; 
u0(N,1:N) = 0.1;  
u0(1:N,1) = 0.1;  
u0(1:N,N) = 0.1;  

%% Solver Loop

% load initial conditions 
u=u0; frame=0; iter=0; t=0; loop=1;

% set time step
dt=cfl*dx^2/k; 

% initilize arrays
Laplace_u=zeros(size(u0));

while t < tEnd
    % time and counter
    t=t+dt; iter=iter+1;
    
    % load old step
    uo=u;
    
    % compute laplace operator
    for i = 2:N-1
        for j = 2:N-1
            Laplace_u(i,j) = (uo(i+1,j)-2*uo(i,j)+uo(i-1,j))/dx^2 + ...
                             (uo(i,j+1)-2*uo(i,j)+uo(i,j-1))/dy^2 ;
        end
    end
    
    % next time step
    u = uo + k*dt*Laplace_u;
    
    % get frame
    if (mod(iter,1)==0) % displays movie frame every 50 time steps
        frame=frame+1;
        surf(x,y,u); view(-45,45);
        h=gca;
        get(h,'FontSize');
        set(h,'FontSize',12);
        xlabel('X','fontSize',12);
        ylabel('Y','fontSize',12);
        colorbar('East');
        title('Heat Diffusion Equation','fontsize',12);
        fh = figure(1);
        set(fh, 'color', 'white');
        F=getframe;
    end
end

%% Postprocessing
% make a movie with the frames
movie(F,frame,1);
