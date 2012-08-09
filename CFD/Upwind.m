%% 1D Linear Advection Equation
%
% du/dt + a du/dx = 0
%
% Subroutine for solving using Upwind Method.
% by Manuel Diaz, manuel.ade'at'gmail.com 
% Institute of Applied Mechanics, 2012.08.08

clear all; close all; clc;

%% Parameters
a = 0.5;
dx = 0.01;
cfl = 0.9;
dt = cfl*dx/a;
t_end = 1.5;

%% Discretization of Domain
x = 1:dx:2;
t = 0:dt:t_end;

%% Initial Condition
n = length(x);
%u_0 = zeros(1,n);
u_0(1:ceil(n/2)) = 2;
u_0((ceil(n/2)+1):n) = 1;

%% Main Loop
u = u_0;
u_next = zeros(1,n);
for k = 1:40
    for j = 2:n-1
        u_next(j) = u(j) - a*dt/dx*(u(j)-u(j-1)); % Upwind
        %u_next(j) = u(j) - dt/(2*dx)*a*(u(j+1)-u(j-1)) ...
        %   +(dt^2/2)/(dx^2)*(a^2)*(u(j+1)-2*u(j)+u(j-1)); % LW
    end
    % BC
    u_next(1) = u_next(2);
    u_next(n) = u_next(n-1);
    %UPDATE info
    u = u_next(1:n);
end
    
%% Plot Results
plot(x,u,'.')
    
