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

%% Simulation Parameters
CFL   = 0.8;    % CFL condition
tEnd  = 0.1;    % out time
theta = 0;      % for FD = -1, MB = 0, BE = 1
quad  = 1;      % for NC = 1 , GH = 2

%% Initialize variables
f     = zeros(ly,lx);
feq   = zeros(ly,lx);
fnext = zeros(ly,lx);

%% Velocity Discretization:
switch quad

    case{1} % Newton Cotes Quadrature:
    V  = [-20,20]; Nv = 200; % 200 cells or 201 points
    vx = linspace(V(1),V(2),Nv+1)';
    dv = vx(2)-vx(1);
    w  = kron(ones(1,Nv/4),[14 32 12 32])'; w(1)=7; w(Nv+1)=7;
        
    case{2} % Gauss Hermite Quadrature:
    Nv = 20;
    [vx,w] = GaussHermite(Nv);
    
    otherwise
        error('Order must be between 1 and 2');
end
    
%% Initial Conditions
% Macroscopic Velocity, Temperature and Fugacity

% Left side:
    u(1,1:lx/2)    = 0;         % velocity
    t(1,1:lx/2)    = 4.383850;  % temperature
    z(1,1:lx/2)    = 0.2253353; % fugacity
% Right side:
    u(1,lx/2+1:lx) = 0;         % velocity
    t(1,lx/2+1:lx) = 8.972544;  % temperature
    z(1,lx/2+1:lx) = 0.1204582; % fugacity

% Equilibrium Distribution function:
%f = 1 /(z^(-1)*exp(()^2/t)+theta);
