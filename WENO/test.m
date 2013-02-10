%% test 
clear all; clc; close all;

% Parameters and ICs
x = linspace(0,1,16);
u = 1*(x<0.5)+0*(x>0.5);
nx = length(x);
dx = max(abs(x(1:end-1)-x(2:end)));
cfl = 0.6;
a = -1; b = 1;
dt = dx*cfl/abs(a);
u1_next = zeros(size(u));
u2_next = zeros(size(u));

%% Test my WENO
h_u = zeros(size(u));
h_d = zeros(size(u));
for i = 3:nx-2
    h_u(i) = WENO_upwind(b*u(i-2:i+2));
    h_d(i) = WENO_downwind(a*u(i-2:i+2));
end


%% Downwind
for j = 3:nx-3
    u1_next(j) = u(j) - dt/dx*(h_d(j+1) - h_d(j));
end

%% Upwind
for j = 4:nx-3
    u2_next(j) = u(j) - dt/dx*(h_u(j) - h_u(j-1));
end
