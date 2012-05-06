%% 2D Heat Transfer D2Q9
% by Manuel Diaz
clc; clear; %close all;

%% Domain Grid
L = 30;    % Length
H = 30;    % Height
n = 30;    % Nodes in x
m = 30;    % Nodes in y
dx = L/n; dy = H/m; dt = 1;
x = 1:dx:L; y = 1:dy:H;
K = 9;      % linkages

%% Initial Condition
f = zeros(m,n,K);
rho = zeros(m,n);    % initial value of the dependent variable rho = T(x,t)
feq = zeros(m,n);

%% Constants
csq   = (dx^2)/(dt^2);
alpha = 0.25;
omega = 1/(3*alpha/(csq*dt)+0.5);
w     = [4/9,1/9,1/9,1/9,1/9,1/36,1/36,1/36,1/36];
cx    = [  0,  1,  0, -1,  0,   1,  -1,  -1,   1];
cy    = [  0,  0,  1,  0, -1,   1,   1,  -1,  -1];
iter  = 200; % time steps

%% Boundary Conditions
twall = 1; 

for i = x;
    for j = y;
        for k = 1:K;
            f(j,i,k) = w(k) * rho(j,i);
            if i == 1;
               f(j,i,k) = w(k) * twall;
            end
        end
    end
end

%% Main Loop
for iterations = 2:iter;
    % Collision process  
    for i=1:K
        rho = rho + f(:,:,k);
    end
        
    for k = 1:K;
        feq = w(k)*rho;
        f(:,:,k) = omega * feq + (1-omega) * f(:,:,k);
    end 
    % Streaming process
    for j = 1:1:m-1;
        for i = 1:1:n-1;
            f(:,j,2) = f( : ,j+1,2);
            f(i,j,6) = f(i+1,j+1,6);
        end
    end
    for j = 1:1:n-1;
        for i = m:1:2;
            f(:,j,3) = f( : ,j+1,3);
            f(i,j,7) = f(i-1,j+1,7);
        end
    end
    for j = n:1:2;
        for i = m:1:2;
            f(:,j,4) = f( : ,j-1,4);
            f(:,j,8) = f(i-1,j-1,8);
        end
    end
    for j = n:1:2;
        for i = 1:1:m-1;
            f(:,j,5) = f( : ,j-1,5);
            f(:,j,9) = f(i+1,j-1,9);
        end
    end
    % Boundary conditions
    f(:,1,2) = w(2)*twall + w(4)*twall - f(:,1,4);
    f(:,1,6) = w(6)*twall + w(7)*twall - f(:,1,7);
    f(:,1,9) = w(9)*twall + w(8)*twall - f(:,1,8);
    
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
    
end

%% Make pretty pictures
contourf(rho)
colormap hot
colorbar('location','southoutside')