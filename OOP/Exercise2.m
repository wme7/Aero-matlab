% Simple Upwind using classes
clear all; clear classes; close all; clc;

% Parameters
a    = 1.5;
cfl  = 0.9;
nx   = 100;
tEnd = 1.0;

% Define structures
domain = Domain;
domain.x     = linspace(0,1,nx);        % xj cell centers
domain.maxdx = abs(max(domain.x(1:nx-1)-domain.x(2:nx)));
domain.un    = zeros(size(domain.x));   % unj cells IC
domain.fu    = a*domain.u;              % flux value @ cell centers

% IC
u0 = sin(2*pi*domain.x);   % u0j cells IC

% time steps
dt = domain.maxdx*cfl/(abs(a));
time = 0:dt:tEnd;

% load Initial condition
domain.u = u0;     % uj cell average values

for tsteps = time
    % Create figs
    range = [0,1,-1,1];
    plot(domain.x,domain.u,'o'); axis(range);
    
    % Upwind Evolution
    dtdx = dt/domain.maxdx;
    ap = 0.5*(a + abs(a));
    an = 0.5*(a - abs(a));
    dup = [0,domain.u(2:nx)] - [0,domain.u(1:nx-1)];
    dun = [domain.u(2:nx),0] - [domain.u(1:nx-1),0];
    domain.un = domain.u - an*dtdx*dun -ap*dtdx*dup ;
    
    %BCs'
    domain.un(end) = domain.un(1);
    domain.un(1) = domain.un(end);
        
    % Update info
    domain.u = domain.un;
            
    % Update figs
    drawnow
end