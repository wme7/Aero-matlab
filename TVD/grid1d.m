function [x,dx] = grid1d(a,b,nx)
%% Description:

% INPUTS
% [a,b]: x range, nx: the desided number of nodes.

% OUTPUTS
% x: vector x's coordinates, dx: cell size.

%% Compute cell size
dx = (b-a)/(nx-1); 

%% Compute vector arrays
x = a:dx:b; 