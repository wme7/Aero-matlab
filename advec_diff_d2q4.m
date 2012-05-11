%% 2D Advection-Diffusion D2Q4
% by Manuel Diaz
clc; clear; %close all;

%% Domain Grid
L = 100;    % Length
H = 100;    % Height
n = 100;    % Nodes in x
m = 100;    % Nodes in y
dx = L/n; dy = H/m; dt = 1;
x = 1:dx:L; y = 1:dy:H;

%% Initial Condition
f1   = zeros(m,n);
f2   = zeros(m,n);
f3   = zeros(m,n);
f4   = zeros(m,n);
feq1 = zeros(m,n);
feq2 = zeros(m,n);
feq3 = zeros(m,n);
feq4 = zeros(m,n);
rho  = zeros(m,n); % initial value of the dependent variable rho = T(x,t)

%% Constants
u     = 0.1;
v     = -0.2; % <- interesting result if we change to -v!
ck    = dx/dt;
csq   = (dx^2)/(dt^2);
alpha = 0.25;
omega = 1/(2*alpha/(dt*csq)+0.5);
w     = [1/4,1/4,1/4,1/4];
cx    = [  1, -1,  0,  0];
cy    = [  0,  0,  1, -1]; 
link  = [  1,  2,  3,  4];
tEnd  = 400; % time steps

%% Movie Parameters
tPlot = 10; frames= tEnd/tPlot; count = 0;
figure(1)
colordef black

%% Boundary Conditions
twall = 1; 

%% Main Loop
for cycle = 2:tEnd;
    % Collision process  
    rho = f1 + f2 + f3 + f4;
    
    feq1 = 0.25 * rho *(1+2*u/ck);
    feq2 = 0.25 * rho *(1-2*u/ck);
    feq3 = 0.25 * rho *(1+2*v/ck);
    feq4 = 0.25 * rho *(1-2*v/ck);
        
    f1 = w(1) * rho;
    f2 = w(2) * rho;
    f3 = w(3) * rho;
    f4 = w(4) * rho;
            
    f1 = omega * feq1 + (1-omega) * f1;
    f2 = omega * feq2 + (1-omega) * f2;
    f3 = omega * feq3 + (1-omega) * f3;
    f4 = omega * feq4 + (1-omega) * f4;
 
    % Streaming process
    f1 = stream2d( f1 , [cx(1),cy(1)]);
    f2 = stream2d( f2 , [cx(2),cy(2)]);
    f3 = stream2d( f3 , [cx(3),cy(3)]);
    f4 = stream2d( f4 , [cx(4),cy(4)]);
    
    % Boundary conditions
    f1(:,1) = twall*w(1) + twall*w(2) - f2(:,1); %Dirichlet BC
    f3(:,1) = twall*w(3) + twall*w(4) - f4(:,1); %Dirichlet BC
    
    f1(m,:) = 0;        %Dirichlet BC
    f2(m,:) = 0;        %Dirichlet BC
    f3(m,:) = 0;        %Dirichlet BC
    f4(m,:) = 0;        %Dirichlet BC
    
    f1(:,n) = 0;        %Dirichlet BC
    f2(:,n) = 0;        %Dirichlet BC
    f3(:,n) = 0;        %Dirichlet BC
    f4(:,n) = 0;        %Dirichlet BC
  
    f1(1,:) = f1(2,:);  %Neumann BC
    f2(1,:) = f2(2,:);  %Neumann BC
    f3(1,:) = f3(2,:);  %Neumann BC
    f4(1,:) = f4(2,:);  %Neumann BC    
    
    % Visualization
    if mod(cycle,tPlot) == 1
        count = count + 1;
        contourf(rho)
        colormap hot
        colorbar('location','southoutside')
        axis tight
        M(count) = getframe;
    end
    
    % Animated gif file
    if mod(cycle,tPlot) == 1
        F = getframe;
        if count == 1
            [im,map] = rgb2ind(F.cdata,256,'nodither');
            im(1,1,1,tEnd/tPlot) = 0;
        end
        im(:,:,1,count) = rgb2ind(F.cdata,map,'nodither');
    end
end

%% Make Movie
movie(M,1,10); % movie(M,n,fps)

%% Export to Gif
imwrite(im,map,'advec_diff_d2q4.gif','DelayTime',0,'LoopCount',3)