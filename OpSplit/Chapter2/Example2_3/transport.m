function u=transport(s,u,t,dx,dir,bndcond)
%% 
% Solves numerically the equation 
% 
% $$u_t + su_x = 0$$
% 
% with different boundary conditions, specified by |bndcond| for a
% time |t|

%% Inital setup
% Some things ensuring that you solve along the columns of u. Must
% transpose if u is a row vector or if dir = 2.
speed=s;
if nargin<6,
	bndcond='Neumann';
end;
if nargin<5,
	dir=1;
elseif (dir>2),
	error('Only 2D transport allowed');
end;
transpose=0;
if size(u,1)==1 || dir==2,
	u=u';
	speed=speed';
	transpose=1;
end;
S=size(u);
speed=reshape(speed,S);
N=S(1); 
%% Solving by a first order  upwind method 
%
% Approximating the solution of u_t + s u_x = 0 by a local upwind
% method. 
%
maxspeed=max(max(abs(speed)));
CFL=0.95;
dt=CFL*dx/maxspeed; 
nstep=ceil(t/dt);
dt=t/nstep;
lambda=dt/dx;
% Setting the inter cell speeds "cmid"
switch lower(bndcond)
	case{'neumann'}
		cext=[speed(1,:);speed;speed(N,:)]; % "Neumann"
	case{'dirichlet'}
		cext=[zeros(1,S(2));speed;zeros(1,S(2))]; % Dirichlet
	case{'periodic'}
		cext=[speed(N,:);speed;speed(1,:)]; % "periodic"
	case {'extrap'}
		cext=[2*speed(1,:)-speed(2,:);speed;2*speed(N,:)-speed(N-1,:)]; % "extrapolated"
end;
cmid=0.5*(cext(1:N+1,:)+cext(2:N+2,:));
ipluss=1.0*(cmid>=0.0); iminus=1.0*(cmid<0.0); 
for j=1:nstep,
	switch lower(bndcond)
	case{'neumann'}
		 uext=[u(1,:);u;u(N,:)];            % "Neumann"
	case{'dirichlet'}
		uext=[zeros(1,S(2));u;zeros(1,S(2))]; % Dirichlet
	case{'periodic'}
		uext=[u(N,:);u;u(1,:)];             % "periodic"
	case {'extrap'}
		uext=[2*u(1,:)-u(2,:);u; 2*u(N,:)-u(N-1,:)];% "extrapolated"
	end;
	umid=uext(1:N+1,:).*ipluss+uext(2:N+2,:).*iminus;
	u=u-lambda*speed.*diff(umid);
end;
if (transpose)
	u=u';
end;
	