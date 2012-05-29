%% BGK-LBM model flow over an obstacle
%
%  XY:            XZ:            YZ:   
%  c8   c2   c7   c12  c5  c11   c16  c5  c15
%    \  |  /        \  |  /        \  |  /  
%  c3 -c19 - c1   c3 -c19 - c1   c3 -c19 - c1
%    /  |  \        /  |  \        /  |  \  
%  c9   c4  c10   c13  c6  c14   c17  c6  c18 
%               
% Using D3Q19 model. At each timestep, particle densities propagate
% outwards in the directions indicated in the figure. An equivalent 
% 'equilibrium' density is found, and the densities relax towards 
% that state, in a proportion governed by omega.(Based on Iain Haslam,
% March 2006 and Jonas Latt 2008.) Modified by Manuel Diaz, May 2012.
% e-mail: manuel.ade@gmail.com

clc; clear; close all;

%% Domain Grid
lx = 100;    % cells in x
ly = 40;    % cells in y
lz = 40;    % cells in z
m  = lx; n = ly; p  = lz; dx = lx/n; dy = ly/m; dt = 1; dz = lz/p;
x  = 1:dx:lx; y = 1:dy:ly; z = 1:dz:lz;
K  = 19;    % Total linkages

%% Obstacles
% Our obstacles are the walls and a single sphere
obst_x = lx/5 + 1;
obst_y = ly/2 + 1;
obst_z = lz/2 + 1;
obst_r = ly/10 + 1;
coly   = [2:(ly-1)];% vector to mark the poinst of inlets and outlets in y
colz   = [2:(lz-1)];% vector to mark the poinst of inlets and outlets in y
in     = 1;         % position of inlet in x
out    = lx;        % position of outlet in x

%% Initial Condition
f = zeros(p,m,n,K);
sum = zeros(p,m,n);
rho = zeros(p,m,n); % initial value of the dependent variable rho = T(x,t)
feq = zeros(p,m,n);

%% Flow Parameters
uMax  = 0.1;
Re    = 100;
nu    = uMax * 2. * obst_r / Re; %kinematic viscority
csq   = (dx^2)/(dt^2);
alpha = 1;
omega = 1./(3*alpha / (csq*dt) + 0.5);

%% Lattice Constants
w     = [4/9,1/9,1/9,1/9,1/9,1/36,1/36,1/36,1/36];
cx    = [ 1 -1  0  0  0  0  1  1  1  1 -1 -1 -1 -1  0  0  0  0  0];
cy    = [ 0  0  1 -1  0  0  1 -1  0  0  1 -1  0  0  1  1 -1 -1  0];
cz    = [ 0  0  0  0  1 -1  0  0  1 -1  0  0  1 -1  1 -1  1 -1  0];
links = [ 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19];
opp   = [ 3  4  1  2  7  8  5  6  7  8 13 14 11 12 17 18 15 16 19];

%% Iteration parameters
tEnd  = 100; % time steps
tPlot = 10;
frames= tEnd/tPlot;
count = 0;

%% Find the cells of our obstacles
[x,y,z] = meshgrid(1:lx,1:ly,1:lz);
obst    = (x - obst_x).^2 + (y - obst_y).^2 + (z - obst_z).^2 <= obst_r^2;
obst([1,ly],:,:) = 4;   % top and bottom bb-cells 
obst(:,:,[1,lz]) = 4;   % left and right bb-cells 
bbRegion = find(obst);  % boolean mask for the bb-cells
% Un comment to view boolean Region:
% PATCH_3Darray(obst,1:ly,1:lx,1:lz)

% L   = y - 2; y_phys = y - 1.5;
% ux  = 4*uMax*(y_phys.*L-y_phys.^2)/L^2;
% uy  = zeros(ly,lx);
% rho = 1;

