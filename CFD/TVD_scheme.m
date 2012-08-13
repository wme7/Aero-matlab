%% 1D Linear Advection Equation
%
% du/dt + a du/dx = 0
%
% Subroutine for solving using TVD Method.
% by Manuel Diaz, manuel.ade'at'gmail.com 
% Institute of Applied Mechanics, 2012.08.12

clear all; close all; clc;

%% Parameters
a = -0.5;
a_m = min(0,a);
a_p = max(0,a);
dx = 0.01;
cfl = 0.9;
dt = cfl*dx/abs(a);
t_end = 0.4;

%% Discretization of Domain
x = 1:dx:2;
t = 0:dt:t_end;

%% Initial Condition
n = length(x);
%u_0 = zeros(1,n);
u_0(1:ceil(n/2)) = 1;
u_0((ceil(n/2)+1):n) = 2;

%% Main Loop

% Initilize vector variables
r = zeros(1,n);
F_rl = zeros(1,n);
F_rh = zeros(1,n);
F_ll = zeros(1,n);
F_lh = zeros(1,n);
F_right = zeros(1,n);
F_left  = zeros(1,n);
u_next = zeros(1,n);

% Load Initial condition
u = u_0;

% Start TVD Loop
for k = t
    for j = 2:n-1
        % smooth measurement factor 'r'
        if u(j) == u(j+1)
            r(j) = 1;
        elseif a > 0
            r(j) = (u(j) - u(j-1)) / (u(j+1) - u(j));
        elseif a < 0
            r(j) = (u(j+2) - u(j+1)) / (u(j+1) - u(j));
        end
        r(1) = 1; r(n) = 1;
    end
        
        % Define Flux Limiter function
        phi = (r + abs(r))./(1 + abs(r));
        
    for j = 2:n-1    
        % Compute fluxes for TVD
        F_rl(j) = a_p*u(j) + a_m*u(j+1);
        F_rh(j) = (1/2)*a*(u(j)+u(j+1)) - (1/2)*(a^2)*(dt/dx)*(u(j+1)-u(j));
        
        F_ll(j) = a_p*u(j-1) + a_m*u(j);
        F_lh(j) = (1/2)*a*(u(j-1)+u(j)) - (1/2)*(a^2)*(dt/dx)*(u(j)-u(j-1));
        
        % Compute next time step
        F_right(j) = F_rl(j) + phi(j)*( F_rh(j) - F_rl(j) );
        F_left(j)  = F_ll(j) + phi(j-1)*( F_lh(j) - F_ll(j) );
        
        u_next(j) = u(j) - dt/dx*(F_right(j) - F_left(j));
    end
    
        % right sided Upwind
        %u_next(j) = u(j) - a*dt/dx*(u(j)-u(j-1)); 
        
        % letf sided Upwind
        %u_next(j) = u(j) - a*dt/dx*(u(j+1)-u(j)); 
        
        % Lax-Wendroff 
        %u_next(j) = u(j) - dt/(2*dx)*a*(u(j+1)-u(j-1)) ...
        %   +(dt^2/2)/(dx^2)*(a^2)*(u(j+1)-2*u(j)+u(j-1)); % LW
        
        % Double sided Upwind
        %u_next(j) = u(j)-a_p*dt/dx*(u(j)-u(j-1))-a_m*dt/dx*(u(j+1)-u(j));
    
        % BC
        u_next(1) = u_next(2);
        u_next(n) = u_next(n-1);
    
        % UPDATE info
        u = u_next(1:n);
end
    
%% Plot Results
plot(x,u,'.'); Ylim([0.5,2.5]); Title('TVD'); 
xlabel('X-Coordinate [-]'); ylabel('U-state [-]');