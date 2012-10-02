%% A finite volume method 
function u=finvol(flux,u,t,dx,dir,bndcond,method,CFL)
%
%  Solves the conservation law u_t + flux(u)_x = 0 by a 
%   finite volume method.
%  "method" is chosen from 'LxF'     : Lax-Friedrichs
%                          'LxW'     : Lax-Wendroff
%                          'McC'     : MacCormac
%                          'upwind'  : upwind 
%
% "flux" is the name of the flux function, must have syntax
%  function f=flux(u,i), where output is the flux if "i" is zero or absent,
%  and the derivative if "i" is one.
%

%% Initial setup
if nargin<8,
	CFL=0.9;
end;
if nargin<7,
	method='LxF';
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
%
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
if strcmp(method,'upwind')
	ms=min(min(speed));
	if ms<0.0,
		display(' Speeds must be positive for upwind method, switching to LxF.');
		method='LxF';
	end;
end;
maxspeed=max(max(abs(speed)));
dt=CFL*dx/maxspeed; 
nstep=ceil(t/dt);
dt=t/nstep;
lambda=dt/dx;
%% The main loop
% ------------------------------------------------------------------
%
% Approximating the solution of u_t + flux(u)_x = 0 by a finite volume method
%

for j=1:nstep,
	switch lower(bndcond)    % Filling in ghost cells according to boundary conditions
	case{'neumann'}
		uext=[u(1,:);u;u(N,:)];                       % "Neumann"
	case{'dirichlet'}
		uext=[zeros(1,S(2));u;zeros(1,S(2))];         % Dirichlet
	case{'periodic'}
		uext=[u(N,:);u;u(1,:)];                       % "periodic"
	case {'extrap'}
		uext=[2*u(1,:)-u(2,:);u; 2*u(N,:)-u(N-1,:)];  % "extrapolated"
	end;
	F=numflux(flux,uext,method,lambda);
	u=u-lambda*diff(F);
end;
if (transpose)
	u=u';
end;
	