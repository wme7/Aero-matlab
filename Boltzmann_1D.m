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
dx = 1/lx;  % delta x
dy = 1/ly;  % delta y

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
u = zeros(1,lx);
t = zeros(1,lx);
z = zeros(1,lx);

% Left side:
    ul = 0;         % velocity
    tl = 4.383850;  % temperature
    zl = 0.2253353; % fugacity
% Right side:
    ur = 0;         % velocity
    tr = 8.972544;  % temperature
    zr = 0.1204582; % fugacity

u(1,1:lx/2) = ul;  u(1,lx/2+1:lx) = ur;
t(1,1:lx/2) = tl;  t(1,lx/2+1:lx) = tr;
z(1,1:lx/2) = zl;  z(1,lx/2+1:lx) = zr;

