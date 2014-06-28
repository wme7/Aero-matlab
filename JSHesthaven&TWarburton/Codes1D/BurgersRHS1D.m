function [rhsu] = BurgersRHS1D(u,epsilon,xL,xR,time)

% function [rhsu] = BurgersRHS1D(u,epsilon,xL, xR, time)
% Purpose  : Evaluate RHS flux in 1D viscous Burgers equation

Globals1D;

% square root of epsilon
Sepsilon = sqrt(epsilon);

% Define field differences at faces
du = zeros(Nfp*Nfaces,K); du(:) = u(vmapM)-u(vmapP);

% impose boundary condition at x=0
%uin =-tanh((xL+0.5-time)/(2*epsilon))+1.0; du(mapI) = 2.0*(u(vmapI)-uin);
%uout= tanh((xR+0.5-time)/(2*epsilon))+1.0; du(mapO) = 2.0*(u(vmapO)-uout);
uin = 1.0;      du(mapI) = 2.0*(u(vmapI)-uin);
uout= 0.1;      du(mapO) = 2.0*(u(vmapO)-uout);

% Compute q and jumps
q = Sepsilon*(rx.*(Dr*u) - LIFT*(Fscale.*(nx.*du/2.0)));
dq = zeros(Nfaces,K); dq(:) = (q(vmapM)-q(vmapP))/2.0;

% impose boundary condition - Dirichlet conditions
dq(mapI) = 0.0; 
dq(mapO) = 0.0;

% Evaluate nonlinear flux
df = zeros(Nfp*Nfaces,K); df(:) = (u(vmapM).^2-u(vmapP).^2)/2.0;

% impose boundary condition
df(mapI)=(u(vmapI).^2- uin.^2); 
df(mapO)=(u(vmapO).^2-uout.^2);

% Compute flux
maxvel = max(max(abs(u)));

% penalty scaling -- See Chapter 7.2
%tau = .25*reshape(N*N./max(2*J(vmapP),2*J(vmapM)), Nfp*Nfaces, K);
tau=0.0;

% flux term
flux = nx.*(df/2.0-Sepsilon*dq) - maxvel/2.0.*du - Sepsilon*tau.*du;

% local derivatives of field
dfdx = rx.*(Dr*(u.^2/2 - Sepsilon*q));

% compute right hand sides of the semi-discrete PDE
rhsu = -(dfdx - LIFT*(Fscale.*flux));
return
