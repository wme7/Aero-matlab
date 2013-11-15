%% Runge Kutta Schemes
%**************************************************************************
% Solve and compare using simple IVP's:
%
%  {1} du/dt = -2u + t + 4, with IC: u(0) = 1
%  {2} du/dt = sin(pi*t),   with IC: u(0) = 1
%  {3} du/dt = EXP(t),      with IC: u(0) = 1
%
% Based on:
%
% [1] Gottlieb and Shu (1998) 
% Mathematics of Computation
% Vol. 67, No. 221. PP. 73-85
%
% [2] Parviz Moin 
% Numerical Method for Engineering and Science
% McGraw-Hill 2nd Ed. (2008)
%
% [3] Jan S. Hesthaven & Tim Warburton
% Nodal Discontinuous Galerkin Methods
% Springer 1st Ed. (2008)
%
% [4] Alex Kanevsky, et al.
% Application of implicit-explicit High order
% Runge-Kutta methods to discontinuous-garlekin methods.
% JCP 225 (2007) 1753-1781
%
% Coded by: Manuel Diaz 2012.12.07
%
%**************************************************************************
%% Clear Workspace
clear all; close all; clc;

%% Define t, dt and tEnd
t0 = 0; tEnd = 2.0; dt = 0.2;
tsteps = t0:dt:tEnd;
t_iter = t0:dt:tEnd-dt;

%% Define number of stages to use
No_stages = 4; % only stages 1, 2, 3 and 4 are available.

%% Choose IVP to study
IVPn = 3;

% Define residual
% Tor this IVP problem the residual is defined as: 
switch IVPn
    case{1}
        residual = @(v,x) -2*v + x + 4;
        range = [t0,tEnd,1,3]; % plot range
    case{2}
        residual = @(v,x) sin(2*pi*x);
        range = [t0,tEnd,1,1.5]; % plot range
    case{3}
        residual = @(v,x) exp(x);
        range = [t0,tEnd,1,9]; % plot range
end

% Compute exact solution
switch IVPn
    case{1}
        u_exact = 1/4*(2*tsteps - 3*exp(-2*tsteps) + 7);
    case{2}
        u_exact = (2*pi+1-cos(2*pi*tsteps))/(2*pi);
    case{3}
        u_exact = exp(tsteps);
end

% Plot exact solution
figure
subplot(1,6,1); plot(tsteps,u_exact); title('Exact'); axis(range);

%% Load IC
u = 1;

%% Run Standard Explicit RK (ERK)
i = 1;
u_ERK(i) = u; % first known value in time
fprintf('ERK\t\t\t');
tic;
for t = t_iter
    switch No_stages
        case{1} % 1st Order Euler time Stepping
            u_next = u + dt*residual(u,t);
            
        case{2} % RK 2nd Order
            u_1 = u + dt/2*residual(u,t);
            u_next = u + dt*residual(u_1,t+1/2*dt);
            
        case{3} % RK 3th Order
            k1 = dt*residual(u,t);
            k2 = dt*residual(u+1/2*k1,t+1/2*dt);
            k3 = dt*residual(u-k1+2*k2,t+dt);
            u_next = u + 1/6*(k1 + 4*k2 + k3);
                        
        case{4} % RK 4th Order
            k1 = dt*residual(u,t);
            k2 = dt*residual(u+1/2*k1,t+1/2*dt);
            k3 = dt*residual(u+1/2*k2,t+1/2*dt);
            k4 = dt*residual(u+k3,t+dt);
            u_next = u + 1/6*(k1 + 2*k2 + 2*k3 + k4);
            
        otherwise
            error('Not a case available')
    end
    i=i+1;
    u_ERK(i) = u_next; % Capture info per iteration
    u = u_next;     % update info
end
toc;

%% Plot ERK solution
subplot(1,6,2); plot(tsteps,u_ERK); title('Standard ERK'); axis(range);

%% Interlude
clear u; u = 1;

%% Run TVD-RK integration
i = 1;
u_TVD(i) = u; % First known value in time
fprintf('TVD-RK\t\t');
tic;
for t = t_iter
    switch No_stages
        case{1} % 1st Order Euler time Stepping
            u_next = u + dt*residual(u,t);
            
        case{2} % TVD-RK 2nd Order
            u_1 = u + dt*residual(u,t);
            u_next = 1/2*u + 1/2*u_1 + 1/2*dt*residual(u_1,t);
            
        case{3} % TVD-RK 3rd Order
            u_1 = u + dt*residual(u,t);
            u_2 = 3/4*u + 1/4*u_1 + 1/4*dt*residual(u_1,t);
            u_next = 1/3*u + 2/3*u_2 + 2/3*dt*residual(u_2,t);
            
        case{4} % TVD-RK 3rd Order
            u_1 = u + dt*residual(u,t);
            u_2 = 3/4*u + 1/4*u_1 + 1/4*dt*residual(u_1,t);
            u_next = 1/3*u + 2/3*u_2 + 2/3*dt*residual(u_2,t);
            
        otherwise
            error('Not a case available')
    end
    i=i+1;
    u_TVD(i) = u_next; % Capture info per iteration
    u = u_next;     % update info
end
toc;

%% Plot TVD-RK solution
subplot(1,6,3); plot(tsteps,u_TVD); title('TVD-RK'); axis(range);

%% Interlude
clear u; u = 1;

%% Run JST-RK 
% general time integrator by Jameson, Schmidt and Turkel:
i = 1;
u_JST(i) = u; %first known value in time
fprintf('JST-RK\t\t');
tic;
for t = t_iter
    %***************************************
    % Stage 0
    sgima = u;
    % Stages 1 to n
    for stage = No_stages:-1:1
       sgima = u + (dt/stage)*residual(sgima,t);
    end
    u_next = sgima;
    %***************************************
    i=i+1;
    u_JST(i) = u_next; % Capture info per iteration
    u = u_next;     % update info
end
toc;

%% Plot JST-RK solution
subplot(1,6,4); plot(tsteps,u_JST); title('RK by JST'); axis(range);

%% Interlude
clear u; u = 1;

%% Run SSP-RK (Strong Stability-Preserving Runge-Kutta Methods) 

i = 1;
u_SSPRK(i) = u; %first known value in time
fprintf('SSP-RK\t\t');
tic;
for t = t_iter
    switch No_stages
        case{1} % 1st Order Euler time Stepping
            u_next = u + dt*residual(u,t);
            
        case{2} % 2nd Order SSP-RK 
            u_1 = u + dt*residual(u,t);
            u_next = 1/2*(u + u_1 + dt*residual(u_1,t+dt));
            
        case{3} % 3rd Order SSP-RK 
            u_1 = u + dt*residual(u,t);
            u_2 = 1/4*(3*u + u_1 + dt*residual(u_1,t+dt));
            u_next = 1/3*(u + 2*u_2 + 2*dt*residual(u_2,t+0.5*dt));
            
        case{4} % 5-stages, 4th-order SSP-RK
            % "It is not possible to construc a fourth-order, four-stage
            % SSP-RK schemes where all coefficients are positive." [3]
            % However, one can derive a fourth-order scheme by allowing a
            % fifth stage. The optimal scheme is given as:
            %
            % Low storage Runge-Kutta coefficients
            rk4a = [            0.0 ...
                    -567301805773.0/1357537059087.0 ...
                    -2404267990393.0/2016746695238.0 ...
                    -3550918686646.0/2091501179385.0  ...
                    -1275806237668.0/842570457699.0];
            rk4b = [ 1432997174477.0/9575080441755.0 ...
                     5161836677717.0/13612068292357.0 ...
                     1720146321549.0/2090206949498.0  ...
                     3134564353537.0/4481467310338.0  ...
                     2277821191437.0/14882151754819.0];
            rk4c = [             0.0  ...
                     1432997174477.0/9575080441755.0 ...
                     2526269341429.0/6820363962896.0 ...
                     2006345519317.0/3224310063776.0 ...
                     2802321613138.0/2924317926251.0];
            resu = 0;
            for s = 1:5
                timelocal = t + rk4c(s)*dt;
                [rhsu] = residual(u,timelocal);
                resu = rk4a(s)*resu + dt*rhsu;
                u = u + rk4b(s)*resu;
            end
            u_next = u;
            
        otherwise
            error('Not a case available')
    end
    i=i+1;
    u_SSPRK(i) = u_next; % Capture info per iteration
    u = u_next;     % update info
end
toc;

%% Plot SSP-RK solution
subplot(1,6,5); plot(tsteps,u_SSPRK); title('SSP-RK'); axis(range);

%% Interlude
clear u; u = 1;

%% Run IMEXRK (Implicit-Explicit Runge-Kutta)
% general time integrator by Jameson, Schmidt and Turkel:
i = 1;
u_ARK(i) = u; %first known value in time
fprintf('ARK\t\t\t');
tic;
 % In progress
toc;
%% Plot TVD-RK solution
%subplot(1,6,6); plot(tsteps,u_ARK); title('ARK'); axis(range);

%% Print Comparisons
error_ERK  = abs(u_exact - u_ERK); fprintf('Standar ERK error\n'); fprintf('%e\t',error_ERK); fprintf('\n\n');
error_TVD = abs(u_exact - u_TVD); fprintf('RK-TVD error\n'); fprintf('%e\t',error_TVD); fprintf('\n\n');
error_JST = abs(u_exact - u_JST); fprintf('JST-TVD error\n'); fprintf('%e\t',error_JST); fprintf('\n\n');
error_SSPRK = abs(u_exact - u_SSPRK); fprintf('JST-SSPRK error\n'); fprintf('%e\t',error_SSPRK); fprintf('\n\n');
%error_ARK = abs(u_exact - u_ARK); fprintf('JST-ARK error\n'); fprintf('%e\t',error_ARK); fprintf('\n\n');

