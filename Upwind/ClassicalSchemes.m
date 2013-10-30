%% Classical Schemes in 1D 
% to solve the scalar advection equation:
%
% du/dt + dF(u)/dx = 0
%
% where $dF/du = a$ i.e. is a constant advection velocity
%
% by Manuel Diaz, manuel.ade'at'gmail.com 
% Institute of Applied Mechanics, 2013.10.12

clear all; close all; clc;

%% Parameters
      a = 0.25;     % Scalar velocity in x direction
     dx = 0.01;     % Spatial step size
    cfl = 0.6;      % Courant Number
     dt = cfl*dx/abs(a); % time step size
   dtdx = dt/dx;    % precomputed to save some flops
  t_end = 0.5;      % End time

%% Discretization of Domain
x = 1:dx:2;     % x grid
n = length(x);  % number of points
t = 0:dt:t_end; % time steps

%% Initial Condition
u0 = ones(1,n);
 x_1 = ceil(2*n/5);
 x_2 = ceil(3*n/5);
u0(x_1:x_2) = 2;

%% 1. Backward Euler
% CAUTION: This method for Hyperbolic conservation laws is unconditionally
% unstable!
u_euler = u0;
% Main loop
for k = t
    for j = 2:n-1
    u_next(j) = u_euler(j) - 1/2*a*dtdx*(u_euler(j+1)-u_euler(j-1));
    end
    u_euler = u_next;
    % BC
    u_euler(1) = u_next(2);
    u_euler(n) = u_next(n-1);
end
u_next = zeros(1,n);

%% 2. Double-Sided Upwind
% Define positive and negative wave velocites
a_p = max(0,a); % a^{+}
a_m = min(0,a); % a^{-}
% load IC
u_upwind = u0;
% Main loop
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

%% 3. Lax-Friedrichs
u_LF = u0;
% Main loop
for k = t
    for j = 2:n-1
    u_next(j) = 1/2*(u_LF(j-1)+u_LF(j+1)) - 1/2*a*dtdx*(u_LF(j+1)-u_LF(j-1));
    end
    u_LF = u_next;
    % BC
    u_LF(1) = u_next(2);
    u_LF(n) = u_next(n-1);
end
u_next = zeros(1,n);

%% 4. Leapfrog
% Leapfrog start in a second time step, therefore we need an initial step,
% we choose backward Euler
u_old = u0;
u_leapfrog = zeros(size(u0));
% Initial stage
for j = 2:n-1
    u_leapfrog(j) = u_old(j) - 1/2*a*dtdx*(u_old(j+1)-u_old(j-1));
end
% BC
u_leapfrog(1) = u_old(2);
u_leapfrog(n) = u_old(n-1);
% leapfrog scheme
for k = t
    for j = 2:n-1
    u_next(j) = u_old(j) - 1/2*a*dtdx*(u_leapfrog(j+1)-u_leapfrog(j-1));
    end
    u_old = u_leapfrog;
    % BC
    u_leapfrog(1) = u_old(2);
    u_leapfrog(n) = u_old(n-1);
    
    u_leapfrog = u_next;
    % BC
    u_leapfrog(1) = u_next(2);
    u_leapfrog(n) = u_next(n-1);
end
u_next = zeros(1,n);

%% 5. Lax-Wendroff
u_LW = u0;
% Main Loop
for k = t
    for j = 2:n-1
    u_next(j) = u_LW(j) - 1/2*a*dtdx*(u_LW(j+1)-u_LW(j-1)) + ...
        1/2*(cfl^2)*(u_LW(j+1)-2*u_LW(j)+u_LW(j-1));
    end
    u_LW = u_next;
    % BC
    u_LW(1) = u_next(2);
    u_LW(n) = u_next(n-1);
end
u_next = zeros(1,n);

%% 6. Beam-Warming
% CAUTION: this is only one sided BW method.
u_BW = u0;
% Main Loop
for k = t
    for j = 3:n
    u_next(j) = u_BW(j) - 1/2*a*dtdx*(3*u_BW(j)-4*u_BW(j-1)+u_BW(j-2)) + ...
        1/2*(cfl^2)*(u_BW(j)-2*u_BW(j-1)+u_BW(j-2));
    end
    u_BW = u_next;
    % BC
    u_BW(2) = u_next(3);
    u_BW(1) = u_next(3);
end
u_next = zeros(1,n);

%% 7. Usign Lax-Wendroff & Beam-Warming
u_combined = u0;
% Main Loop
for k = t
    for j = 3:n-1
    u_next(j) = u_combined(j) ...
        - 1/4*a*dtdx*(3*u_combined(j)-4*u_combined(j-1)+u_combined(j-2)) ...
        + 1/4*(cfl^2)*(u_combined(j)-2*u_combined(j-1)+u_combined(j-2)) ...
        - 1/4*a*dtdx*(u_combined(j+1)-u_combined(j-1)) ...
        + 1/4*(cfl^2)*(u_combined(j+1)-2*u_combined(j)+u_combined(j-1));
    end
    u_combined = u_next;
    % BC
    u_combined(2) = u_next(3);
    u_combined(1) = u_next(3);
    u_combined(n) = u_next(n-1);
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
subplot(241);
hold on
plot(x,u_euler,'.'); plot(x,u_exact,'k'); ylabel('u(x,t)'); ylim([0.5,2.5]); title 'Backward Euler'; 
hold off
subplot(242);
hold on
plot(x,u_upwind,'.'); plot(x,u_exact,'k'); ylabel('u(x,t)'); ylim([0.5,2.5]); title 'Double-Sided Upwind'; 
hold off
subplot(243);
hold on
plot(x,u_LF,'.'); plot(x,u_exact,'k'); ylabel('u(x,t)'); ylim([0.5,2.5]); title 'Lax-Friedrichs';
hold off
subplot(244);
hold on
plot(x,u_leapfrog,'.'); plot(x,u_exact,'k'); ylabel('u(x,t)'); ylim([0.5,2.5]); title 'Leapfrog';
hold off
subplot(245);
hold on
plot(x,u_LW,'.'); plot(x,u_exact,'k'); ylabel('u(x,t)'); ylim([0.5,2.5]); title 'Lax-Wendroff';
hold off
subplot(246);
hold on
plot(x,u_BW,'.'); plot(x,u_exact,'k'); ylabel('u(x,t)'); ylim([0.5,2.5]); title 'Beam-Warming'; xlabel('x');
hold off
subplot(247);
hold on
plot(x,u_combined,'.'); plot(x,u_exact,'k'); ylabel('u(x,t)'); ylim([0.5,2.5]); title 'Combine: LW & BW'; xlabel('x');
hold off
