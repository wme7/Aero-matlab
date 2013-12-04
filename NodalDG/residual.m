function dF = residual(u,flux,dflux,Lift,Dr)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Compute the Residual for 1d wave equation using Nodal DG 
%
%                       Residual = dF/dxi 
%                 where F = is our Flux function
%
%              coded by Manuel Diaz, NTU, 2013.10.29
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% compute fluxes in node coordinates
f = flux(u); 

% Interpolate u and flux values at the boundaries of Ij
% use only: {'LGL','ChebyshevMod'}
u_lbd = u(1,:);
u_rbd = u(end,:);
f_lbd = f(1,:);
f_rbd = f(end,:);

% Build Numerical fluxes across faces
u_pface = [u_lbd,0]; % + side 
u_nface = [0,u_rbd]; % - side 

% Apply Periodic BCs
u_nface(1) = u_nface(end); % left BD
u_pface(end) = u_pface(1); % right BD

% Apply Neumann BCs
%u_nface(1) = u_pface(1); % left BD
%u_pface(end) = u_nface(end);% u_pface(end); % right BD

% LF numerical flux
alpha = max(max(abs(dflux(u)))); 
nflux = 0.5*(flux(u_nface)+flux(u_pface)-alpha*(u_pface-u_nface));
nfluxL = nflux(1:end-1); nfluxR = nflux(2:end);

% Compute the derivate: F = f + lL*(nfluxL-f_bdL) - lR*(nfluxR-f_bdR)
dF = -Dr*f + Lift*[(nfluxL-f_lbd);(f_rbd-nfluxR)];

dF = -dF;
