function u_0 = u_zero2d(x,y)
%% Initial Condition (IC)
% This subroutine creates an special IC for testing purposes.
% The present IC models a rectangular domain with four equally sized
% regions with diferrent initial state value, u, an adimensional value.
%
%  Region 1 = 1.0 [-]
%  Region 2 = 0.0 [-]
%  Region 3 = 0.5 [-]
%  Region 4 = 1.5 [-]
%
% by Manuel Diaz, manuel.ade'at'gmail.com 
% Institute of Applied Mechanics, 2012.09.06
%% Initial condition
% Parameters
nx = length(x);
ny = length(y);

% Preallocate u_0,
u_0 = zeros(ny,nx);

% Parameters of regions dimensions
x_middle = ceil(nx/2);
y_middle = ceil(ny/2);
l_1 = 1:x_middle; l_2 = x_middle+1:nx;
h_1 = 1:y_middle; h_2 = y_middle+1:ny;

% Initial Condition for our 2D domain
u_0(h_1,l_1) = 0.70; % region 1
u_0(h_1,l_2) = 0.10; % region 2
u_0(h_2,l_1) = 0.90; % region 3
u_0(h_2,l_2) = 0.50; % region 4