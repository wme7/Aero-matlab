%% Discontinuous Galerkin (DG) 1D subroutine
% to solve the scalar advection equation:
%
% du/dt + dF(u)/dx = 0
%
% where $F(u) = au$ or $dF/du = a$
%
% by Manuel Diaz, manuel.ade'at'gmail.com 
% Institute of Applied Mechanics, 2012.09.03

clear all; close all; clc; Globals1D;

%% Routine main parameters
      a = 0.5;    % advection speed
     Nv = 4;     % Number of elements en V_h space
      k = 5;      % order of polynomial used for approximation
    CFL = 0.6;    % Courant Number
  t_end = 0.5;    % Final time
NODETOL = 1e-10;  % Node Tolerance
     Np = k+1;    % Number of points per element
    Nfp = 1;      % Number of faces per point
 Nfaces = 2;      % Number of faces
    
%% Discretization of Domain
% Grid Parameters:
      v = [-1,1];         % Domain range
     dv = (v(2)-v(1))/Nv; % Spatial step size
     dt = CFL*dv/abs(a);  % Time step size
   dtdx = dt/dv;          % precomputed to save some flops
 Nsteps = ceil(t_end/dt); % Number of time steps
% Element Grid Generation:
     vx = v(1):dv:v(2);   % Elements grid points
      t = 0: dt :t_end;   % time steps
     Nx = Nv+1;           % Total number of elements grid points
% Provide info of connectivity between elements:
 EToV = zeros(Nv,2);
for i = 1:Nv;  EToV(i,1) = i;  EToV(i,2) = i+1; end 

%% Compute basic mathematical objects
% Compute Legendre-Gauss-Lobatto nodes and weigths, and 
r = JacobiGL(0,0,k);

% Compute Vandermonde matrix (P) and it's inverse (invP):
% V_{ij} = phi_j(r_i)
V  = Vandermonde1D(k, r); invV = inv(V);

% Compute Dr matrix
Dr = Dmatrix1D(k, r, V);

%% LIFT
% Create surface integral terms
Emat = zeros(Np,Nfaces*Nfp);

% Define Emat
Emat(1,1) = 1.0; Emat(Np,2) = 1.0;

% inv(mass matrix)*\s_n (L_i,L_j)_{edge_n}
LIFT = V*(V'*Emat);

%% X-coordinate and maps
% build coordinates of all the nodes
va = EToV(:,1); vb = EToV(:,2);
 x = ones(k+1,1)*vx(va) + 0.5*(r+1)*(vx(vb)-vx(va));

% calculate geometric factors
[rx,J] = GeometricFactors1D(x,Dr);

% Compute masks for edge nodes
fmask1 = find( abs(r+1) < NODETOL)'; 
fmask2 = find( abs(r-1) < NODETOL)';
Fmask  = [fmask1;fmask2]';
Fx = x(Fmask(:), :);

% Build surface normals and inverse metric at surface
[nx] = Normals1D();
Fscale = 1./(J(Fmask,:));

% Build connectivity matrix
[EToE, EToF] = Connect1D(EToV);

% Build connectivity maps
[vmapM, vmapP, vmapB, mapB] = BuildMaps1D;
 
%% IC
% Initial Condition
u_0 = zeros(1,Nv*Np);
% jump = [0.9 1.1];
% for i = 1:Nv*Np
%     if jump(1) <= x(i) && x(i) <= jump(2)
%         u_0(i) = 2;
%     else
%         u_0(i) = 1;
%     end
% end
% Gaussian Jump
u0 = 0.1 + 0.5*exp(-20*x.^2);
% Load Initial Condition
u = reshape(u_0,Np,Nv);
%u = sin(x);

%% Main Loop
% Runge-Kutta residual storage  
resu = zeros(Np,Nv); 

time = 0;
% outer time step loop 
for tstep=1:Nsteps
    for INTRK = 1:5
        timelocal = time + rk4c(INTRK)*dt;
        [rhsu] = AdvecRHS1D(u, timelocal, a);
        resu = rk4a(INTRK)*resu + dt*rhsu;
        u = u+rk4b(INTRK)*resu;
    end;
    % Increment time
    time = time+dt;
end;

%% Compute exact Solution
% Displacement = a*time;
u_exact = zeros(1,Nv*Np);
jump = [0.9+a*time,1.1+a*time];
for i = 1:Nv*Np
    if jump(1) <= x(i) && x(i) <= jump(2)
        u_exact(i) = 2;
    else
        u_exact(i) = 1;
    end
end
u_exact = reshape(u_exact,Np,Nv);