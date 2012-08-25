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
a = 0.5; 
a_p = max(0,a); a_m = min(0,a);
dx = 0.01;
cfl = 0.8;
dt = cfl*dx/abs(a);
dtdx = dt/dx; % precomputed to save some flops
t_end = 0.2;
limiter = 2;

%% Discretization of Domain
x = 1:dx:2;
t = 0:dt:t_end;

%% Initial Condition
n = length(x);
u_0 = ones(1,n);
 x_1 = ceil(2*n/5);
 x_2 = ceil(3*n/5);
u_0(x_1:x_2) = 2;

%% Main Loop
% Load initial condition into the domain
u = u_0; 

% Initialize vector variables 
u_next = zeros(1,n);
theta  = zeros(1,n);
r = zeros(1,n);
F_rl = zeros(1,n);
F_rh = zeros(1,n);
F_ll = zeros(1,n);
F_lh = zeros(1,n);
F_right = zeros(1,n);
F_left  = zeros(1,n);

for k = t
    % Compute the smoothness, r(j), from the data, u(j).
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
    
    % Compute the Flux Limiter
    switch limiter
        case{1} % Van Leer
            phi = (r + abs(r))./(1 + abs(r));
        case{2} % Superbee
            %phi = max(0,min(2.*r,0),min(r,2));
            phi_star = max(min(2.*r,0),min(r,2));
            phi = max(0,phi_star);
        case{3} % Minmod
            phi = max(0,min(1,r));
        case{4} % koren
            %phi = max(0,min(2*r,(2/3)*r+1/3,2));
            phi_star = min(2*r,(2/3)*r+1/3);
            phi_hat  = min(phi_star,2);
            phi = max(0,phi_hat);
        otherwise 
            error('Limiters available: 1, 2, 3, and 4')
    end
    
    for j = 2:n-1    
        % Compute fluxes for TVD
        F_rl(j) = a_p*u(j) + a_m*u(j+1);
        F_rh(j) = (1/2)*a*(u(j)+u(j+1)) - (1/2)*(a^2)*dtdx*(u(j+1)-u(j));
        F_right(j) = F_rl(j) + phi(j)*( F_rh(j) - F_rl(j) );
        
        F_ll(j) = a_p*u(j-1) + a_m*u(j);
        F_lh(j) = (1/2)*a*(u(j-1)+u(j)) - (1/2)*(a^2)*dtdx*(u(j)-u(j-1));
        F_left(j)  = F_ll(j) + phi(j-1)*( F_lh(j) - F_ll(j) );
        
        % Compute next time step
        u_next(j) = u(j) - dtdx*(F_right(j) - F_left(j));
    end
    
    % BC
    u_next(1) = u_next(2);
    u_next(n) = u_next(n-1);
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
    u_next(j) = u_LW(j) - 1/2*cfl*(u_LW(j+1)-u_LW(j-1)) + ...
        1/2*(cfl^2)*(u_LW(j+1)-2*u_LW(j)+u_LW(j-1));
    end
    u_LW = u_next;
    % BC
    u_LW(1) = u_next(2);
    u_LW(n) = u_next(n-1);
end
u_next = zeros(1,n);

%% Exact Solution
x_exact_1 = 1.4 + a*t_end;
x_exact_2 = 1.6 + a*t_end;
u_exact = zeros(1,n);
for j = 1:n
    if x(j) >= x_exact_1 && x(j) <= x_exact_2
    u_exact(j) = 2;
    else
    u_exact(j) = 1;
    end
end
    
%% Plot Results
subplot(311);
hold on
plot(x,u_upwind,'.'); plot(x,u_exact,'k'); xlabel('X-Coordinate [-]'); ylabel('U-state [-]'); yLim([0.5,2.5]); title 'Double-Sided Upwind'; 
hold off
subplot(312);
hold on
plot(x,u_LW,'.'); plot(x,u_exact,'k'); xlabel('X-Coordinate [-]'); ylabel('U-state [-]'); yLim([0.5,2.5]); title 'Lax-wendroff';
hold off
subplot(313);
hold on
plot(x,u,'.'); plot(x,u_exact,'k'); xlabel('X-Coordinate [-]'); ylabel('U-state [-]'); yLim([0.5,2.5]); title 'TVD';
hold off

