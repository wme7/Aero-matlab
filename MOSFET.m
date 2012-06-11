%% Escalar Poisson Equation 
% 
% Routine for solving a PNP-MOSFET junction for the region:
% ________________________
%|++++|    ---N---   |++++|
%|+P++|              |+P++|
%|++++|              |++++|
%|____|              |____|
%|                        |
%|                        |
%|                        |
%|________________________|

%% Governing equation
%
% $$\frac{\partial}{\partial x}\left (k(x,y)\frac{\partial{\phi }}{\partial{x}}\right)+
%   \frac{\partial}{\partial y}\left (k(x,y)\frac{\partial{\phi}}{\partial{y}}\right)=
%   \tau (\hat{n}-\hat{n}_d)$$
%
% Where: 
%
% $\hat{n}$ is the number density of a doped region at time t,
% $\hat{n}_d$ is the initial number density of the doped regions and $\tau$
% is a constant. 

clear, clc, %close all;

%% IC's
omega = 1.4;
delta = 0.1;

L  = 1.8; dx = delta; x = 1:dx:L; n  = 18;
H  = 1.2; dy = delta; y = 1:dy:H; m  = 12;

tau = 1;

n_1 = 1; % initial doping for region 1 (silicon)
n_2 = 3; % initial doping for region 2 (gate A)
n_3 = 3; % initial doping for region 3 (gate B)

%% Grid
n_d = ones(m,n);            % Silicon
n_d(1:m/4, 1:n/6 ) = 3;     % gate A
n_d(1:m/4, 5*n/6+1:n) = 3;  % gate B
phi = zeros(m,n);           % electric potential grid
phi_tilde = zeros(m,n); phi_next = zeros(m,n);

% x gates locations:
x_g1 = 1:n/6;
x_gc = 2*n/6+1:4*n/6;
x_g2 = 5*n/6+1:n;

% x silicon location at the surface 
x_s1 = 1*n/6+1:2*n/6;
x_s2 = 4*n/6+1:5*n/6;

%% Computing k(x,y) on the grid:
k = ones(m,n); % Assuming all k values constants

%% BC's
% bolean mask for gates locations 
x_gates = zeros(1,n); 
x_gates(x_g1) = 1; x_gates(x_gc) = 1; x_gates(x_g2) = 1;

% Dirichlet:
phi(1,x_g1) =  0.0;     % gate A
phi(1,x_gc) = -0.8;     % central gate
phi(1,x_g2) =  1.0;     % gate B
% Neumann BC: Every where else!

%% Main Loop
A = zeros(m,n); B = zeros(m,n); C = zeros(m,n); D = zeros(m,n);

% we know the boundary values from the begining
phi_tilde(1,x_g1) =  0.0; % gate A
phi_tilde(1,x_gc) = -0.8; % central gate
phi_tilde(1,x_g2) =  1.0; % gate B

% r = 0.0001; iter = 0; r_iter = 1; 
% while r_iter >= r
for s=1:1100
    for i=2:m-1
        for j=2:n-1
            % Update the change in k varaible
            A(i,j)=k(i,j)+k(i+1,j);
            B(i,j)=k(i,j)+k(i,j+1);
            C(i,j)=k(i,j)+k(i-1,j);
            D(i,j)=k(i,j)+k(i,j-1);
            if j == 2 % Matrix's left column, Neumann BC
                phi_tilde(i,j)=...
                    1/(k(i+1,j)+k(i-1,j)+3*k(i,j)+k(i,j+1))*(...
                    A(i,j)*phi(i+1,j)+...
                    B(i,j)*phi(i,j+1)+...
                    C(i,j)*phi_tilde(i-1,j));
            elseif j == n-1 % Matrix's right column, Neumann BC
                phi_tilde(i,j)=...
                    1/(k(i+1,j)+k(i-1,j)+3*k(i,j)+k(i,j-1))*(...
                    A(i,j)*phi(i+1,j)+...
                    C(i,j)*phi_tilde(i-1,j)+...
                    D(i,j)*phi_tilde(i,j-1));
            elseif i == 2 && x_gates(j)==1
                % Matrix Upper row, for the regions with Neumann BC
                phi_tilde(i,j)=...
                    1/(k(i+1,j)+3*k(i,j)+k(i,j+1)+k(i,j-1))*(...
                    A(i,j)*phi(i+1,j)+...
                    B(i,j)*phi(i,j+1)+...
                    D(i,j)*phi_tilde(i,j-1));
            elseif i == m-1 % Matrix Lower row, Neumann BC
                phi_tilde(i,j)=...
                    1/(k(i-1,j)+3*k(i,j)+k(i,j+1)+k(i,j-1))*(...
                    B(i,j)*phi(i,j+1)+...
                    C(i,j)*phi_tilde(i-1,j)+...
                    D(i,j)*phi_tilde(i,j-1));
            else
                phi_tilde(i,j)=...
                    1/(k(i+1,j)+k(i-1,j)+4*k(i,j)+k(i,j+1)+k(i,j-1))*(...
                    A(i,j)*phi(i+1,j)+...
                    B(i,j)*phi(i,j+1)+...
                    C(i,j)*phi_tilde(i-1,j)+...
                    D(i,j)*phi_tilde(i,j-1));
                    %+(delta^2)/4*;
            end
        end
    end
    phi_tilde(2:m-1,1) = phi_tilde(2:m-1, 2 ); % Left matrix column, Neumann BC
    phi_tilde(2:m-1,n) = phi_tilde(2:m-1,n-1); % Right matrix column, Neumann BC
    phi_tilde(1, x_s1) = phi_tilde(2, x_s1  ); % Upper regions matrix row regions, Neumann BC
    phi_tilde(1, x_s2) = phi_tilde(2, x_s2  );
    phi_tilde(m,2:n-1) = phi_tilde(m-1,2:n-1); % Lower matrix column, Neumann BC
    % OSR steps:
    phi_next = phi + omega*( phi_tilde - phi );
    r_iter = norm(phi_next-phi,2);
    phi = phi_next;
    %iter = iter+1;
end
%fprintf('iterations: %4.1f \n',iter);

%% Plot Figures
figure
contourf(phi)
title(['SOR Method, h =',num2str(delta),', \omega =',num2str(omega)])
%title(['SOR Method, h =',num2str(h),', \omega =',num2str(omega),', iterations: ',num2str(iter)])
xlabel('x points')
ylabel('y points')
colormap hot
colorbar('location','southoutside')
