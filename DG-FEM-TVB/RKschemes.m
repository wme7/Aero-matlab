%% Runge Kutta Schemes
%**************************************************************************
% Solve and compare two simple IVP's:
%
%  {1} du/dt = -2u + t + 4, with IC: u(0) = 1
%  {2} du/dt = sin(pi*t),   with IC: u(0) = 1
%  {3} du/dt = EXP(t),      with IC: u(0) = 1
%
% Based on:
%
% Gottlieb and Shu (1998) 
% Mathematics of Computation
% Vol. 67, No. 221. PP. 73-85
%
% Parviz Moin 
% Numerical Method for Engineering and Science
% McGraw-Hill 2nd Ed. (2008)
%
% Coded by: Manuel Diaz 2012.12.07
%
%**************************************************************************
%% Clear Workspace
clear all; close all; clc;

%% Define t, dt and tEnd
t0 = 0; tEnd = 1.0; dt = 0.2;
tsteps = t0:dt:tEnd;
t_iter = t0:dt:tEnd-dt;

%% Define number of stages to use
No_stages = 4; % only stages 1,2 and 3 are available.

%% Choose IVP to study
IVPn = 3;

% Define residual
% Tor this IVP problem the residual is defined as: 
switch IVPn
    case{1}
        residual = @(v,x) -2*v + x + 4;
    case{2}
        residual = @(v,x) sin(2*pi*x);
    case{3}
        residual = @(v,x) exp(x);
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
range = [t0,tEnd,1,3]; % plot range 
subplot(1,4,1); plot(tsteps,u_exact); title('Exact'); axis(range);

%% Load IC
u = 1;

%% Run RK Classical integration schemes
i = 1;
u_RK(i) = u; % first known value in time
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
    u_RK(i) = u_next; % Capture info per iteration
    u = u_next;     % update info
end
%% Plot TVD-RK solution
subplot(1,4,2); plot(tsteps,u_RK); title('Classic RK'); axis(range);

%% Interlude
clear u; u = 1;

%% Run TVD-RK integration
i = 1;
u_TVD(i) = u; %first known value in time
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

%% Plot TVD-RK solution
subplot(1,4,3); plot(tsteps,u_TVD); title('TVD-RK'); axis(range);

%% Interlude
clear u; u = 1;

%% Run JST-RK 
% general time integrator by Jameson, Schmidt and Turkel:
i = 1;
u_JST(i) = u; %first known value in time
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

%% Plot TVD-RK solution
subplot(1,4,4); plot(tsteps,u_JST); title('RK by JST'); axis(range);

%% Print Comparisons
error_RK  = abs(u_exact - u_RK); fprintf('Classic RK error\n'); fprintf('%e\t',error_RK); fprintf('\n\n');
error_TVD = abs(u_exact - u_TVD); fprintf('RK-TVD error\n'); fprintf('%e\t',error_TVD); fprintf('\n\n');
error_JST = abs(u_exact - u_JST); fprintf('JST-TVD error\n'); fprintf('%e\t',error_JST); fprintf('\n\n');


