function dF = residual(u,L,dg,flux,dflux,quad)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Compute the Residual for 1d wave equation using CPR/FR 
%
%                       Residual = dF/dxi 
%             where F = is our Correcting Flux function
%
%              coded by Manuel Diaz, NTU, 2013.10.29
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% compute fluxes in node coordinates
f = flux(u);

% Interpolate u and flux values at the boundaries of Ij
switch quad
    case 'LGL'
        u_lbd = u(1,:);
        u_rbd = u(end,:);
        f_lbd = f(1,:);
        f_rbd = f(end,:);
    otherwise
        u_lbd = L.lcoef*u;
        u_rbd = L.rcoef*u;
        f_lbd = L.lcoef*f;
        f_rbd = L.rcoef*f;
end
% Build Numerical fluxes acroos faces
u_pface = [u_lbd,0]; % + side
u_nface = [0,u_rbd]; % - side

% Apply Periodic BCs
u_nface(1) = u_nface(end); % left BD
u_pface(end) = u_pface(1); % right BD

% LF numerical flux
alpha = max(max(abs(dflux(u))));
nflux = 0.5*(flux(u_nface)+flux(u_pface)-alpha*(u_pface-u_nface));
nfluxL = nflux(1:end-1); nfluxR = nflux(2:end);

% flux derivate
df = L.dcoef*f;

% Compute the derivate: F = f + gL*(nfluxL-f_bdL) + gR*(nfluxR-f_bdR)
dF = df + dg.RR*(nfluxL - f_lbd) + dg.RL*(nfluxR - f_rbd);

