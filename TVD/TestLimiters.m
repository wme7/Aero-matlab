%% Total Variation Diminishing (TVD) 1D Test subroutine
% to solve the scalar advection equation:
%
% du/dt + dF(u)/dx = 0
%
% where $dF/du = a$
%
% Using flux limiter functions: 
% Cases: {1} Van Leer
%        {2} Superbee
%        {3} Minmod
%        {4} Koren
%
% by Manuel Diaz, manuel.ade'at'gmail.com 
% Institute of Applied Mechanics, 2012.09.06

clear all; close all; clc;

%% Parameters
      a = -0.5;      % Scalar velocity in x direction
    a_p = max(0,a); % a^{+}
    a_m = min(0,a); % a^{-}
     dx = 0.01;     % Spatial step size
    cfl = 0.9;      % Courant Number
     dt = cfl*dx/abs(a); % time step size
   dtdx = dt/dx;    % precomputed to save some flops
  t_end = 0.2;      % End time
%limiter = 1;        % Options: 1(Vl), 2(Sb), 3(Mm), 4(koren)

%% Discretization of Domain
x = 1:dx:2;     % x grid
n = length(x);  % number of points
t = 0:dt:t_end; % time steps

%% IC
u_0 = u_zero1d(x);

%% 1. TVD with Van Leer limiter 
limiter = 1;
% Load initial condition into the domain
u_vl = u_0; 
for k = t
    % Compute the smoothness factors, r(j), from data, u(j).
    [r] = theta1d(u_vl,a);
    
    % Compute the Flux Limiter
    [phi] = fluxlimiter1d(r,limiter);
    
    % Compute TVD Fluxes
    [F_left,F_right] = TVDflux1d(u_vl,a,dtdx,phi);
        
    % Compute next time step
    u_next = u_vl - dtdx*(F_right - F_left);
        
    % BC
    u_next(1) = u_next(2);      %Neumann BC
    u_next(n) = u_next(n-1);  	%Neumann BC
    
    % UPDATE info
    u_vl = u_next(1:n);
end
clear r phi F_left F_right u_next;

%% 2. TVD with Superbee limiter 
limiter = 2;
% Load initial condition into the domain
u_sp = u_0; 
for k = t
    % Compute the smoothness factors, r(j), from data, u(j).
    [r] = theta1d(u_sp,a);
    
    % Compute the Flux Limiter
    [phi] = fluxlimiter1d(r,limiter);
    
    % Compute TVD Fluxes
    [F_left,F_right] = TVDflux1d(u_sp,a,dtdx,phi);
        
    % Compute next time step
    u_next = u_sp - dtdx*(F_right - F_left);
        
    % BC
    u_next(1) = u_next(2);      %Neumann BC
    u_next(n) = u_next(n-1);    %Neumann BC
    
    % UPDATE info
    u_sp = u_next(1:n);
end
clear r phi F_left F_right u_next;


%% 3. TVD with minmod limiter 
limiter = 3;
% Load initial condition into the domain
u_mm = u_0; 
for k = t
    % Compute the smoothness factors, r(j), from data, u(j).
    [r] = theta1d(u_mm,a);
    
    % Compute the Flux Limiter
    [phi] = fluxlimiter1d(r,limiter);
    
    % Compute TVD Fluxes
    [F_left,F_right] = TVDflux1d(u_mm,a,dtdx,phi);
        
    % Compute next time step
    u_next = u_mm - dtdx*(F_right - F_left);
        
    % BC
    u_next(1) = u_next(2);      %Neumann BC
    u_next(n) = u_next(n-1);    %Neumann BC
    
    % UPDATE info
    u_mm = u_next(1:n);
end
clear r phi F_left F_right u_next;

%% Exact Solution
u_exact = u_exact1d(x,a,t_end,dx);

%% Plot Pretty figures
subplot(311);
hold on
plot(x,u_vl,'.'); plot(x,u_exact,'k'); xlabel('X-Coordinate [-]'); ylabel('U-state [-]'); yLim([0.5,2.5]); title 'Van Leer'; 
hold off
subplot(312);
hold on
plot(x,u_sp,'.'); plot(x,u_exact,'k'); xlabel('X-Coordinate [-]'); ylabel('U-state [-]'); yLim([0.5,2.5]); title 'Superbee';
hold off
subplot(313);
hold on
plot(x,u_mm,'.'); plot(x,u_exact,'k'); xlabel('X-Coordinate [-]'); ylabel('U-state [-]'); yLim([0.5,2.5]); title 'Minmod';
hold off