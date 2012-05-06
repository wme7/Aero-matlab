%% 2D Heat Transfer D2Q4
% by Manuel Diaz
clc; clear; %close all;

%% Domain Grid
L = 30;    % Length
H = 30;    % Height
n = 30;    % Nodes in x
m = 30;    % Nodes in y
dx = L/n; dy = H/m; dt = 1;
x = 1:dx:L; y = 1:dy:H;

%% Initial Condition
f1 = zeros(m,n);
f2 = zeros(m,n);
f3 = zeros(m,n);
f4 = zeros(m,n);
rho= zeros(m,n); % initial value of the dependent variable rho = T(x,t)
feq= zeros(m,n);

%% Constants
csq   = (dx^2)/(dt^2);
alpha = 0.25;
omega = 1/(2*alpha/(dt*csq)+0.5);
w     = [1/4,1/4,1/4,1/4];
iter  = 200; % time steps

%% Boundary Conditions
twall = 1; 

%% Main Loop
for k = 2:iter;
    % Collision process  
    rho = f1 + f2 + f3 + f4;
    % for this case k1=k2=k3=k4=1/4 then feq1=...feq4=feq
    feq = 0.25 * rho; 
        
    f1 = w(1) * rho;
    f2 = w(2) * rho;
    f3 = w(3) * rho;
    f4 = w(4) * rho;
            
    f1 = omega * feq + (1-omega) * f1;
    f2 = omega * feq + (1-omega) * f2;
    f3 = omega * feq + (1-omega) * f3;
    f4 = omega * feq + (1-omega) * f4;
    % Streaming process
    for i = 1:1:n-1;
        f1(i,:) = f1(i+1,:);  % Streaming
    end
    for i = n:-1:2;
        f2(i,:) = f2(i-1,:);  % Streaming        
    end
    for j = 1:1:m-1;
        f3(:,j) = f3(:,j+1);  % Streaming
    end
    for j = m:-1:2;
        f4(:,j) = f4(:,j-1);  % Streaming
    end
    % Boundary conditions
    f1(:,1) = twall*w(1) - f2(:,1);    %Dirichlet BC
    f3(:,1) = twall*w(2) - f4(:,1);    %Dirichlet BC
    
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
end

%% Make pretty pictures
contourf(rho)
colormap hot
colorbar('location','southoutside')