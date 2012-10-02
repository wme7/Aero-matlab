function u=LxF(flux,u,t,dx,dir,bndcond)
%% Lax-Friedrichs mehod 
% for a scalar conservation law.
% Some things ensuring that you solve along the columns of u. Must
% transpose if u is a row vector or if dir = 2.
%
% "flux" is the name of the flux function, must have syntax
%  
%  function f=flux(u,i), where output is the flux if "i" is zero or absent,
%  and the derivative if "i" is one.
%

%% Initial setup
if nargin<6,
	bndcond='Neumann';
end;
if nargin<5,
	dir=1;
end;
if (dir>2),
	error('You can easily extend this to 3D yourself ...');
end;
transpose=0;
if (size(u,1)==1) || (dir==2),
	u=u';
	transpose=1;
end;
S=size(u);
N=S(1); 
speed=feval(flux,u,1);
maxspeed=max(max(abs(speed)));
CFL=0.475;
dt=CFL*dx/maxspeed; 
nstep=ceil(t/dt);
dt=t/nstep;
lambda=dt/dx;

%% Iteration loop

for j=1:nstep,
    %  start by filling in "artificial" cells at the end of the strip
    %  according to boundary condition.
	switch lower(bndcond)
	case{'neumann'}
		 uext=[u(1,:);u;u(N,:)];                     % "Neumann"
	case{'dirichlet'}
		uext=[zeros(1,S(2));u;zeros(1,S(2))];       % Dirichlet
	case{'periodic'}
		uext=[u(N,:);u;u(1,:)];                     % "periodic"
	case {'extrap'}
		uext=[2*u(1,:)-u(2,:);u; 2*u(N,:)-u(N-1,:)];% "extrapolated"
	end;
	fext=feval(flux,uext);
    % LxF evaluation
	u=0.5*(uext(3:N+2,:)+uext(1:N,:)-lambda*(fext(3:N+2,:)-fext(1:N,:)));
end;
if (transpose)
	u=u';
end;
	