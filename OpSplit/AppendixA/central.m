%% A central scheme 
function u=central(flux,u,t,dx,dir,bndcond,CFL,limiter)
%
%  Solves the conservation law u_t + flux(u)_x = 0 by a Nessayu-Tadmor
%  central  method
%  using the limiter function "limiter"
%
% "flux" is the name of the flux function, must have syntax
%  function f=flux(u,i), where output is the flux if "i" is zero or absent,
%  and the derivative if "i" is one.
%

%% Initial setup
if nargin<8,
  limiter='minmod';
end;
if nargin<7,
  CFL=0.5;
end;
if nargin<6,
  bndcond='Neumann';
end;
if nargin<5,
  dir=1;
end;
if (dir>2),
  error('You can easily extend this to 3D yourself ...');
end;

% Some things ensuring that you solve along the columns of u. Must
% transpose if u is a row vector or if dir = 2.
%
transpose=0;
S=size(u);
N=S(1);
if (N==1||dir==2),
  u=u';
  transpose=1;
end;
S=size(u);
N=S(1); 
speed=feval(flux,u,1);
maxspeed=max(max(abs(speed)));
dt=CFL*dx/maxspeed; 
nstep=ceil(t/dt);
if (mod(nstep,2)) % 'Correcting' the staggered grid...
    nstep=nstep+1;
end;
dt=t/nstep;
lambda=dt/dx; 
ll=dx/(8*dt);
%% The main loop 
for j=1:nstep,
  switch lower(bndcond)   % Inserting ghost cells
   case{'neumann'}
    uext=[u(1,:);u;u(N,:)];                       % "Neumann"
	uuext=[uext(1,:);uext;uext(N+2,:)];       
   case{'dirichlet'}
    uext=[zeros(1,S(2));u;zeros(1,S(2))];         % Dirichlet
	uuext=[zeros(1,S(2));uext;zeros(1,S(2))];
   case{'periodic'}
    uext=[u(N,:);u;u(1,:)];                       % "periodic"
	uuext=[uext(N+2,:);uext;uext(1,:)];        
   case {'extrap'}
    uext=[2*u(1,:)-u(2,:);u; 2*u(N,:)-u(N-1,:)];  % "extrapolated"
	uuext=[2*uext(1,:)-uext(2,:);u; 2*uext(N+2,:)-uext(N+1,:)]; 
  end;
  ff=feval(flux,uuext);
  dff=diff(ff);
  sigma=feval(limiter,dff(1:N+2,:),dff(2:N+3,:));
  duu=diff(uuext);
  s=feval(limiter,duu(1:N+2,:),duu(2:N+3,:));
  uh=uext-0.5*lambda*sigma;
  g=feval(flux,uh)+ll*s;
  uh=0.5*(uext(1:N+1,:)+uext(2:N+2,:))-lambda*diff(g);
  if mod(j,2),
	  u=uh(2:N+1,:);
  else
	 u=uh(1:N,:);
  end;
end;
% u=0.5*(uh(1:N,:)+uh(2:N+1,:)); % Uncomment this is you want some
%smoothing... 
if (transpose)
  u=u';
end;
	