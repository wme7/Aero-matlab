function [x,dx,y,dy] = grid2d(a,b,nx,c,d,ny)
%% Description:

% INPUTS
% [a,b]: x range, nx: the desided number of nodes.
% [c,d]: y range, ny: the desided number of nodes.

% OUTPUTS
% x: vector x's coordinates, dx: cell size.
% y: vector y's coordinates, dy: cell size.

%% Compute cell size
dx = (b-a)/(nx-1); dy = (d-c)/(ny-1);

%% Compute vector arrays
x = a:dx:b; y = c:dy:d; 