%% Total Variation Diminishing (TVD) 1D subroutine
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
% Institute of Applied Mechanics, 2012.08.12

clear all; close all; clc;

%% Parameters
      a = 0.5;      % Scalar velocity in x direction
    a_p = max(0,a); % a^{+}
    a_m = min(0,a); % a^{-}
     dx = 0.01;     % Spatial step size
    cfl = 0.8;      % Courant Number
     dt = cfl*dx/abs(a); % time step size
   dtdx = dt/dx;    % precomputed to save some flops
  t_end = 0.1;      % End time
limiter = 1;        % Options: 1(Vl), 2(Sb), 3(Mm), 4(koren)

%% Discretization of Domain
x = 1:dx:2;     % x grid
n = length(x);  % number of points
t = 0:dt:t_end; % time steps

%% Initial Condition
u_0 = ones(1,n);
 x_1 = ceil(2*n/5);
 x_2 = ceil(3*n/5);
u_0(x_1:x_2) = 2;

%% Main Loop
% Load initial condition into the domain
u = u_0; 

% Initialize vector variables 
u_next = zeros(1,n);

for k = t
    % Compute the smoothness factors, r(j), from data, u(j).
    [r] = theta1d(u,a);
    
    % Compute the Flux Limiter
    [phi] = fluxlimiter1d(r,limiter);
    
    % Compute TVD Fluxes
    [F_left,F_right] = TVDflux1d(u,a,dtdx,phi);
        
    % Compute next time step
    u_next = u - dtdx*(F_right - F_left);
        
    % BC
    u_next(1) = u_next(2);   % left boundary condition (fixed flow value)
    u_next(n) = u_next(n-1); % right boundary condition (fixed flow value)
    
    % UPDATE info
    u = u_next(1:n);
end

%% 2. Double-Sided Upwind
u_upwind = u_0;
for k = t
    for j = 2:n-1
    u_next(j) = u_upwind(j) - a_p*dt/dx*(u_upwind(j) - u_upwind(j-1)) ...
        - a_m*dt/dx*(u_upwind(j+1) - u_upwind(j));
    end
    u_upwind = u_next;
    % BC
    u_upwind(1) = u_next(2);
    u_upwind(n) = u_next(n-1);
end
u_next = zeros(1,n);

%% 3. Lax-Wendroff
u_LW = u_0;
for k = t
    for j = 2:n-1
    u_next(j) = u_LW(j) - 1/2*a*dt/dx*(u_LW(j+1)-u_LW(j-1)) + ...
        1/2*(cfl^2)*(u_LW(j+1)-2*u_LW(j)+u_LW(j-1));
    end
    u_LW = u_next;
    % BC
    u_LW(1) = u_next(2);
    u_LW(n) = u_next(n-1);
end
u_next = zeros(1,n);

%% Exact Solution
% Displacement:
x_add = ceil(a*t_end/dx);
% Computing exact solution:
u_exact = ones(1,n);
 x_1 = ceil(2*n/5) + x_add;
 x_2 = ceil(3*n/5) + x_add;
u_exact(x_1:x_2) = 2;
    
%% Plot Results
subplot(311);
hold on
plot(x,u_upwind,'.'); plot(x,u_exact,'k'); xlabel('X-Coordinate [-]'); ylabel('U-state [-]'); ylim([0.5,2.5]); title 'Double-Sided Upwind'; 
hold off
subplot(312);
hold on
plot(x,u_LW,'.'); plot(x,u_exact,'k'); xlabel('X-Coordinate [-]'); ylabel('U-state [-]'); ylim([0.5,2.5]); title 'Lax-wendroff';
hold off
subplot(313);
hold on
plot(x,u,'.'); plot(x,u_exact,'k'); xlabel('X-Coordinate [-]'); ylabel('U-state [-]'); ylim([0.5,2.5]); title 'TVD';
hold off