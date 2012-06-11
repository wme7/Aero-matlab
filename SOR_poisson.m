%% Escalar Poisson Equation rutine solver for
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
% $$\frac{\partial^2 \phi}{\partial x^2}+\frac{\partial^2 phi}{\partial x^2}=
%   -\tau (\hat{n}-\hat{n}_d)$$
%
% Where: 
%
% $\hat{n}$ is the number density of a doped region at time t,
% $\hat{n}_d$ is the initial number density of the doped regions and $\tau$
% is a constant. 

clear, clc, %close all;

%% IC's
omega = 1.4;    % Relaxation constant
delta = 0.1;    % using an uniform grid in x and y, h = delta
tau = 1;        % any constant

L = 18; dx = delta; x = 1:dx:L; n  = round(L/delta);
H = 12; dy = delta; y = 1:dy:H; m  = round(H/delta);

%% Grid
n_1 = 1; % initial doping value for silicon
n_2 = 3; % initial doping value for gate 1
n_3 = 3; % initial doping value for gate 2

n_i = zeros(m,n);           % array for silicon doping chaging in every iteration
n_d = ones(m,n);            % initial silicon doping array
n_d(1:m/4, 1:n/6 )   = n_2; % set values of n_d for gate 1
n_d(1:m/4, 5*n/6+1:n)= n_3; % set values of n_d for gate 2

phi = zeros(m,n);           % Electric potential grid (our goal)
phi_tilde = zeros(m,n); phi_next = zeros(m,n); % Arrays for SOR

% gates locations in x:
x_g1 = 1:n/6;
x_gc = 2*n/6+1:4*n/6;
x_g2 = 5*n/6+1:n;

% silicon location at the surface in x:
x_s1 = 1*n/6+1:2*n/6;
x_s2 = 4*n/6+1:5*n/6;

%% BC's
% bolean mask for gates locations 
x_gates = zeros(1,n); 
x_gates(x_g1) = 1; x_gates(x_gc) = 1; x_gates(x_g2) = 1;

% Dirichlet:
phi(1,x_g1) =  0.0;     % gate A
phi(1,x_gc) = -0.8;     % central gate
phi(1,x_g2) =  1.0;     % gate B
% Neumann BC: Every where else!

% we know this boundary values of next 'tilde' step from the begining
phi_tilde(1,x_g1) =  0.0; % gate A
phi_tilde(1,x_gc) = -0.8; % central gate
phi_tilde(1,x_g2) =  1.0; % gate B

%% Main Loop

%r = 0.0001; iter = 0; r_iter = 1; 
%while r_iter >= r
%n_i = n_d; % evaluate laplace formulation!

for s = 1000 % for fixed number of iterations
    % Using a sqrt func to increase n_i until it reaches 2
    if s < 10
        n_i = n_d.*(2-1/s^2); 
        %n_i = n_d.*(2-1/iter^2); 
    else
        n_i = 2.*n_d;
    end
    % SOR: Computing phi_tilde loop:
    for i=2:m-1
        for j=2:n-1
            if j == 2 % Matrix's left column, Neumann BC
                phi_tilde(i,j) = 1/3*(...
                    phi(i+1,j)+...
                    phi(i,j+1)+...
                    phi_tilde(i-1,j))...
                    +1/3*(delta^2)*tau*(n_i(i,j)-n_d(i,j));
            elseif j == n-1 % Matrix's right column, Neumann BC
                phi_tilde(i,j) = 1/3*(...
                    phi(i+1,j)+...
                    phi_tilde(i-1,j)+...
                    phi_tilde(i,j-1))...
                    +1/3*(delta^2)*tau*(n_i(i,j)-n_d(i,j));
            elseif i == 2 && x_gates(j)==1
                % Matrix Upper row, for the regions with Neumann BC
                phi_tilde(i,j) = 1/3*(...
                    phi(i+1,j)+...
                    phi(i,j+1)+...
                    phi_tilde(i,j-1))...
                    +1/3*(delta^2)*tau*(n_i(i,j)-n_d(i,j));
            elseif i == m-1 % Matrix Lower row, Neumann BC
                phi_tilde(i,j) = 1/3*(...
                    phi(i,j+1)+...
                    phi_tilde(i-1,j)+...
                    phi_tilde(i,j-1))...
                    +1/3*(delta^2)*tau*(n_i(i,j)-n_d(i,j));
            else
                phi_tilde(i,j) = 1/4*(...
                    phi(i+1,j)+...
                    phi(i,j+1)+...
                    phi_tilde(i-1,j)+...
                    phi_tilde(i,j-1))...
                    +1/4*(delta^2)*tau*(n_i(i,j)-n_d(i,j));
            end
        end
    end
    % Update for Neumann BC's
    phi_tilde(2:m , 1) = phi_tilde(2:m , 2  ); % Left array column
    phi_tilde(2:m , n) = phi_tilde(2:m , n-1); % Right array column
    phi_tilde(1, x_s1) = phi_tilde(2, x_s1  ); % Upper array regions 
    phi_tilde(1, x_s2) = phi_tilde(2, x_s2  );
    phi_tilde(m,2:n-1) = phi_tilde(m-1,2:n-1); % Lower array column
    
    % SOR: compute phi next using realation constant: omega
    phi_next = phi + omega*( phi_tilde - phi );
    %r_iter = norm(phi_next-phi);
    phi = phi_next;
    %iter = iter+1;
end
%fprintf('iterations: %4.1f \n',iter);

%% Plot Figures
contourf(phi)
title(['SOR Method, h = ',num2str(delta),', \omega = ',num2str(omega),', iterations: ',num2str(s)])
xlabel('x points')
ylabel('y points')
colormap hot
colorbar('location','southoutside')