%% 2D Heat Transfer D2Q9
% by Manuel Diaz 05/2012
clc; clear; close all;

%% Domain Grid
L = 10;    % Length
H = 10;    % Height
n = 10;    % Nodes in x
m = 10;    % Nodes in y
dx = L/n; dy = H/m; dt = 1;
x = 1:dx:L; y = 1:dy:H;
K = 9;     % linkages

%% Initial Condition
f = zeros(m,n,K);
sum = zeros(m,n);
rho = zeros(m,n);    % initial value of the dependent variable rho = T(x,t)
feq = zeros(m,n);

%% Constants
csq   = (dx^2)/(dt^2);
alpha = 0.25;
omega = 1/(3*alpha / (csq*dt) + 0.5);
w     = [4/9,1/9,1/9,1/9,1/9,1/36,1/36,1/36,1/36];
cx    = [  0,  1,  0, -1,  0,   1,  -1,  -1,   1];
cy    = [  0,  0,  1,  0, -1,   1,   1,  -1,  -1];
links = [  1,  2,  3,  4,  5,   6,   7,   8,   9];
tEnd  = 100; % time steps
tPlot = 10;
frames= tEnd/tPlot;
count = 0;

%% Boundary Conditions
twall = 1; 

for k = links;
    
    if k == 1; % the central node
        %f(:,:,k) = w(k) * twall;
    else
        f(:,:,k) = w(k) * rho;
    end
end

%% Main Loop
for cycle = 1:tEnd;
    % Collision process
    for k = links;
        sum = sum + f(:,:,k);
    end
    rho = sum;
    for k = links;
        feq = w(k)*rho;
        f(:,:,k) = omega*feq + (1-omega)*f(:,:,k);
    end 
    
    % Streaming process
    for k = links
        f(:,:,k) = stream2d(f(:,:,k) , [cx(k),cy(k)]);
    end
    
    % Boundary conditions
    f(:,1,2) = w(2)*twall + w(4)*twall - f(:,1,4);
    f(:,1,6) = w(6)*twall + w(8)*twall - f(:,1,8);
    f(:,1,9) = w(9)*twall + w(7)*twall - f(:,1,7);
    
    f(:,n,4) = -f(:,n,2);
    f(:,n,7) = -f(:,n,9);
    f(:,n,8) = -f(:,n,6); 
    
    f(m,:,8) = -f(m,:,6);
    f(m,:,5) = -f(m,:,3);
    f(m,:,9) = -f(m,:,7);

    f(1,:,1) = f(2,:,1);
    f(1,:,2) = f(2,:,2);
    f(1,:,3) = f(2,:,3);
    f(1,:,4) = f(2,:,4);
    f(1,:,5) = f(2,:,5);
    f(1,:,6) = f(2,:,6);
    f(1,:,7) = f(2,:,7);
    f(1,:,8) = f(2,:,8);
    f(1,:,9) = f(2,:,9);
    
    % Visualization
    if mod(cycle,tPlot) == 1
        count = count + 1;
        contourf(rho)
        colormap hot
        colorbar('location','southoutside')
        title '2D Heat_Equation using LBM D2Q9'
        M(count)=getframe;
    end
end

%% Make Movie
movie(M,2,10); % movie(M,n,fps)