%% 1D Boltzmann Equation
% Numerical solution of the boltzmann equation to recover Euler macroscopic
% continuum solution. By Manuel Diaz 2012.10.6
%
% Boltzmann Trasnport Equation:
%
% $$\frac{\partial f}{\partial t}+\vec F\cdot \nabla_p f + \vec v
% \cdot\nabla_{\vec x} f =\widehat{\Omega } (f)$$
%
clc;  clear all;  close all;

%% Simulation Parameters
CFL    = 0.8;    % CFL condition
tEnd   = 0.1;    % Out time
theta  = 0;      % for FD = -1, MB = 0, BE = +1
quad   = 1;      % for NC = 1 , GH = 2
method = 1;      % for TVD = 1, WENO3 = 2, WENO5 = 3
input  = 5;      % Reimann IC case

%% Space Discretization
X  = linspace(0,1,8);          % Physical domain -x
Y  = linspace(0,2,8);          % Physical domain -y
dx = max(X(2:end)-X(1:end-1)); % delta x
dy = max(Y(2:end)-Y(1:end-1)); % delta y
[x,y] = meshgrid(X,Y);

%% Velocity Discretization:
switch quad

    case{1} % Newton Cotes Quadrature:
    V  = [-10,10];  % range: a to b
    lv = 40;        % nodes desired (may not the actual value)
    [v,w,k] = cotes_xw(V(1),V(2),lv,5); % cotes Degree 5
        
    case{2} % Gauss Hermite Quadrature:
    lv = 20;        % nodes desired (the actual value)
    [v,w] = GaussHermite(lv); % for integrating range: -inf to inf
     k = 1;         % quadrature constant.
    
    otherwise
        error('Order must be between 1 and 2');
end
%% Velocity-Space Grid:
% The actual lv and lx value will be computed using 'lenght' vector function:
clear lv; lv = length(v);
clear lx; lx = length(x);
[nv,nx] = meshgrid(1:lv,1:lx); % get coordinate of matrix indices
% Initialize Arrays
u = zeros(lv,lx);  t = zeros(lv,lx);   r=zeros(lv,lx);    
v = v*ones(1,lx);  w = w*ones(1,lx);
f = zeros(lv,lx);  feq = zeros(lv,lx); fnext = zeros(lv,lx);

%% Initial Conditions
% Load Macroscopic Velocity, Temperature and Fugacity
[r,u,w,t] = reimann_IC2d(x,y,input);

%% Distribution function:
% Compute initial Equilibrium distribution function:
    feq = 1 ./ ( exp( (v-u).^2 ./ t) ./ r  + theta );
% Plot initial Equilibrium function:
    %surface(feq)

%% Recover Initial Macroscopic Properties:
% Using Quadrature rules to integrate for:
    n   = k*sum(w .* feq);    % Number density
    u_n = k*sum(v .* w .* feq);   % Macrospic velocity in x
    E   = k*sum(1/2*( v.^2 ).* w .* feq); % Energy
% Plot Macroscopic Properties of feq:
    %plot(x,n,x,u,x,E) 

%% Advection Scheme
% By negleting any force field acting over our domian, The Boltzmann equation
% would resemble to a pure advection equation. Thus we use WENO or TVD
% subroutine to compute the advection of the information in the grid:
switch advec
    
    case{1} %TVD 0(h^2)
        
    case{2} %WENO k = 3 i.e. O(h^5) 
        
    otherwise
        error('Order must be between 1 and 2');
end

%% Recover Macroscopic Properties:
% Using Quadrature rules to integrate for:
    n   = k*sum(w .* feq);    % Number density
    u_n = k*sum(v .* w .* feq);   % Macrospic velocity in x
    E   = k*sum(1/2*( v.^2 ).* w .* feq); % Energy
% Plot Macroscopic Properties of feq:
     plot(x,n,x,u,x,E) 