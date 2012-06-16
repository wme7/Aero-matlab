%% 1D Boltzmann Equation
% Numerical solution of the boltzmann equation neglecting any force acting
% over the fisical domain. 
clc;  clear all;  close all;

%% Position Grid
L  = 1;     % Physical domain length
H  = 1;     % Physical domain height
lx = 12;    % number of cells in x direction
ly = 12;    % number of cells in y direction
[y,x] = meshgrid(1:ly,1:lx); % get coordinate of matrix indices
dx = L/lx;  % delta x
dy = H/ly;  % delta y

%% Velocity Grid
V  = [-20,20];
dv = (V(1)-V(2))/20;
v  = V(1):dv:V(2);

%% Simulation Parameters
CFL  = 0.8;     % CFL condition
tEnd = 0.1;     % out time

%% Initialize variables
f   = zeros(ly,lx);
feq = zeros(ly,lx);
fnext = zeros(ly,lx);

%% Gauss Hermite Polynomials
[nx,w] = GaussHermite(20);

%% Initial Conditions
% Left side:
    u(1,1:lx/2)    = 0;         % velocity
    t(1,1:lx/2)    = 4.383850;  % temperature
    z(1,1:lx/2)    = 0.2253353; % fugacity
% Right side:
    u(1,lx/2+1:lx) = 0;         % velocity
    t(1,lx/2+1:lx) = 8.972544;  % temperature
    z(1,lx/2+1:lx) = 0.1204582; % fugacity


