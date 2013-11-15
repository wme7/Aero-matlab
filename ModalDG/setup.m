function [x,xi,w,V] = setup(k,x_nodes,quadn)
% Inputs:
% k: Element polynomial degree that will be use in every cell/element
% x: Row vector with cell's nodes. i.e. ( x_{j} , x_{j+1} , ... )
% quadn: strategy to build information

% Define grids center points. i.e. ( x_{j+1/2} , x_{j+3/2} , ... )
x_mid = x_nodes(1:end-1)+(x_nodes(2:end)-x_nodes(1:end-1))/2; % cells middle points

%% Parameters
nx = length(x_mid); % Number of elements (and/or center nodes)
bb = ones(1,nx);    % row boolean operator
dx = abs(x_nodes(2:end)-x_nodes(1:end-1)); % individual cell's size 
np = k+1;           % number for points per element

% Build cells/elements using the central points
x = repmat(x_mid,k+1,1); % Create cells

%% Main Routine
% Build elements points, weigths and compute the vandermonde matrix
switch quadn
    case{1} % Scaled Guass-Lobatto quadrature weights and abscissas 
        [xi,w]   = sGaussLobatto(np); % scaled to [-1/2, 1/2]
        xj       = xi*dx;
        w        = w*bb;
        x        = x + xj; % Create points in every cells
        V        = sLegMat(k,xi);
        
    case{2} % Guass-Lobatto quadrature weights and abscissas 
        [xi,w,V] = GaussLobatto(np); % for the interval [-1,1]
        xj       = xi*dx/2;
        w        = w*bb;
        x        = x + xj; % Create points in every cells
        
    case{3} % Guass-Legendre quadrature weights and abscissas
        [xi,w,V] = GaussLegendre(np); % for the interval [-1,1]
        xj       = xi*dx/2;
        w        = w*bb;
        x        = x + xj; % Create points in every cells
                
    case{4} % Guass-Legendre quadrature weights and abscissas 
        [xi,w,V] = GaussRadau(np); % for the interval [-1,1]
        xj       = xi*dx/2;
        w        = w*bb;
        x        = x + xj; % Create points in every cells
        
    otherwise
        error('quadrature not available')
end