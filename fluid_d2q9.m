%% 2D Lattice Boltzmann (BGK) model for ideal fluid.
%
%  c7  c3   c6  Using D2Q9 model. At each timestep, particle densities propagate
%    \  |  /    outwards in the directions indicated in the figure. An
%  c4 -c1 - c2  equivalent 'equilibrium' density is found, and the densities
%    /  |  \    relax towards that state, in a proportion governed by omega.
%  c8  c5   c9  Based on Iain Haslam, March 2006 and Jonas Latt 2008.
%               Modified by Manuel Diaz, May 2012.
%               e-mail: manuel.ade@gmail.com

clc; clear; close all;

%% Domain Grid
lx = 20;  % cells in x
ly = 10;  % cells in y
m  = lx; n = ly; dx = lx/n; dy = ly/m; dt = 1;
x = 1:dx:lx; y = 1:dy:ly;
K = 9;     % total linkages

%% Obstacles
% Our obstacles are the walls and a single sphere
obst_x = lx/5 + 1;
obst_y = ly/2 + 1;
obst_r = ly/10 + 1;
col    = [2:(ly-1)];% vector to mark the poinst of inlets and outlets in y
in     = 1;         % position of inlet in x
out    = lx;        % position of outlet in x

%% Initial Condition
f = zeros(m,n,K);
sum = zeros(m,n);
rho = zeros(m,n);    % initial value of the dependent variable rho = T(x,t)
feq = zeros(m,n);

%% Flow Parameters
uMax  = 0.1;
Re    = 100;
nu    = uMax * 2. * obst_r / Re; %kinematic viscority
csq   = (dx^2)/(dt^2);
alpha = 1;
omega = 1./(3*alpha / (csq*dt) + 0.5);

%% Lattice Constants
w     = [4/9,1/9,1/9,1/9,1/9,1/36,1/36,1/36,1/36];
cx    = [  0,  1,  0, -1,  0,   1,  -1,  -1,   1];
cy    = [  0,  0,  1,  0, -1,   1,   1,  -1,  -1];
links = [  1,  2,  3,  4,  5,   6,   7,   8,   9];
opp   = [  1,  4,  5,  2,  3,   8,   9,   6,   7];

%% Iteration parameters
tEnd  = 100; % time steps
tPlot = 10;
frames= tEnd/tPlot;
count = 0;

%% Find the cells of our obstacles
[x,y]    = meshgrid(1:lx,1:ly);
obst     = (x - obst_x).^2 + (y - obst_y).^2 <= obst_r^2;
obst([1,ly],:) = 1;     % top and bottom bb-cells 
bbRegion = find(obst);  % boolean mask for the bb-cells

L   = y - 2; y_phys = y - 1.5;
ux  = 4*uMax*(y_phys.*L-y_phys.^2)/L^2;
uy  = zeros(ly,lx);
rho = 1;

