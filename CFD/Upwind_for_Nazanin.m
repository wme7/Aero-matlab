%% 1D Linear Advection Equation
% for Nazanin
%
% du/dt + a du/dx = 0
%
% Subroutine for solving using Upwind Method.
% by Manuel Diaz, manuel.ade'at'gmail.com 
% Institute of Applied Mechanics, 2013.03.08

clear all; close all; clc;

%% Parameters
dx  = 1/200; % using 200 points
cfl  = 0.9;
t_end = 0.5;

%% Discretization of Domain
x = 0:dx:1;
n = length(x);

%% Initial Condition
u0 = zeros(1,n);
idx = find( x>=0.2 & x<=0.4);
u0(idx) = 1;

%% Main Loop
% Initilize vector variables
u_next = zeros(size(u0)); % u^{n+1}

% Load Initial condition
u = u0;

% Main Loop
t = 0;
while t < t_end
    % plot every time step
    plot(x,u,'.')
    
    % Compute advection speed
    a = (1+x.^2)./(1+2*x*t+2*x.^2+x.^4);
    a_m = min(0,a);
    a_p = max(0,a);
    
    % compute time step
    dt = cfl*dx/abs(max(a));
    t = t + dt;
    
    for j = 2:n-1    
        % Using double Sided Upwind
        u_next(j) = u(j) - a_p(j)*dt/dx*(u(j)-u(j-1)) - a_m(j)*dt/dx*(u(j+1)-u(j));
    end
    % BC
    u_next(1) = 0;              % Dirichlet BC
    %u_next(n) = u_next(n-1);    % Neumann BC
    
    % UPDATE info
    u = u_next;
    
    % update figure
    pause(0.01); drawnow
end

%% Compute exact solution at t = t_end
u_exact = zeros(1,n);
idx = find( x>=(0.2 - t_end/(1+0.2^2) ) & x<=(0.4 - t_end/(1+0.2^2)));
u_exact(idx) = 1;

%% Plot final Result
hold on
plot(x,u,'.')
%plot(x,u_exact,'-k')
hold off
