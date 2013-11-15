function [rhsu] = iBurgersRHS1D(u)

% function [rhsrho, rhsrhou, rhsEner] = EulerRHS1D(rho, rhou ,Ener)
% Purpose  : Evaluate RHS flux in 1D Euler

Globals1D;

% compute maximum velocity for LF flux
lm = abs(u);

% Compute fluxes
uf = u.^2/2; 

% compute source term
us = u.^2;

% Compute jumps at internal faces
du  =zeros(Nfp*Nfaces,K);  du(:)  =  u(vmapM) -  u(vmapP); 
duf =zeros(Nfp*Nfaces,K); duf(:)  = uf(vmapM) - uf(vmapP);
LFc =zeros(Nfp*Nfaces,K); LFc(:)  = max(lm(vmapP),lm(vmapM));

% Compute fluxes at interfaces
duf(:) = nx(:).*duf(:)/2.0-LFc(:)/2.0.*du(:); 

% Boundary conditions 

% Dirichlet
% uin    = 0.0;   
% uout   = 0.0;   
% 
% % Set fluxes at inflow/outflow
% ufin =uin.^2/2; 
% lmI=lm(vmapI)/2; nxI=nx(mapI);
% duf (mapI)=nxI*(uf (vmapI)-ufin )/2.0-lmI*(u(vmapI) -uin);  
% 
% ufout=uout.^2/2; 
% lmO=lm(vmapO)/2; nxO=nx(mapO);
% duf (mapO)=nxO*(uf(vmapO) - ufout)/2.0-lmO*(u(vmapO)- uout);  

% Neumann
duf (mapI) = 0;
duf (mapO) = 0;

% compute right hand sides of the PDE's
rhsu  = -rx.*(Dr*uf)  + LIFT*(Fscale.*duf) + ((MassMatrix^-1)*us);
return
