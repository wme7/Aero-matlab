%% Scalar Advection with WENO
clear all; clc; close all;

% Parameters and ICs
x = linspace(0,1,16);
u = 1*(x<0.5)+0*(x>0.5);
nx = length(x);
dx = max(abs(x(1:end-1)-x(2:end)));
cfl = 0.6;
a = -1;
dt = dx*cfl/abs(a);

% Initialize Arrays
%u1_next = zeros(size(u));
%u2_next = zeros(size(u));
u_next = zeros(size(u));

%% Flux Spliting
fluxsplit = 1;
switch fluxsplit
    case{1} % Godunov
        f = a*u;
        fp = 0.5*(f + abs(f));
        fn = 0.5*(f - abs(f));
    case{2} % Using Lax friedrichs
        f = a*u; alpha = max(abs(a));
        fp = 0.5*(f + alpha*u); %flux^{+}
        fn = 0.5*(f - alpha*u); %flux^{-}
    otherwise 
        error('only cases 1 and 2 are available')
end

%% Test my WENO
h_u = zeros(size(u));
h_d = zeros(size(u));
for i = 3:nx-2
    h_u(i) = WENO_upwind(fp(i-2:i+2));
    %h_d(i) = WENO_downwind(fn(i-2:i+2)); 
    h_d(i-1) = WENO_downwind(fn(i-2:i+2)); 
    % why i-1?: Draw the full domain to see why!
end

hh = [h_u;h_d];
h = sum(hh);

%% Downwind
% for j = 3:nx-4
%     u1_next(j) = u(j) - dt/dx*(h_d(j+1) - h_d(j));
% end 

%% Upwind
% for j = 4:nx-3
%     u2_next(j) = u(j) - dt/dx*(h_u(j) - h_u(j-1));
% end

%% Putting both together!
for j = 4:nx-3
    u_next(j) = u(j) - dt/dx*(h(j) - h(j-1));
end
