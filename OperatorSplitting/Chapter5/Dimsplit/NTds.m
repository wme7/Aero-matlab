function U=NTds(U,X,Y,xflux,yflux,N,T,varargin)
%%
%   Performs dimensional splitting on the 2D scalar conservation law
%   U_t + xflux(U)_x + yflux(U)_y = 0 with initial data U.
%   Uses Strang splitting with timestep dt=T/N. The domain is the square
%   spanned by X and Y. varargin may contain 'wbar' if waitbar is desired,
%   'periodic' followd by [xmin xmax], [ymin ymax]. The one-dimensional
%   equations are solved using a second order central scheme. 
%

%%
% Initial setup
addpath ../../AppendixA;
wb=0; 
nv=length(varargin);
i=1;
boundary_cond='neumann';
while i<=nv,
	if strcmp(varargin(i),'wbar'),
		wb=1;
		i=i+1;
	elseif strcmp(varargin(i),'periodic'),
		boundary_cond='periodic';
		i=i+1;
	else
      display('Unknown option in NTds');
      i=i+1;
	end;
end;
dt=T/N;
xdir=2; ydir=1;
dx=X(2)-X(1);
dy=Y(2)-Y(1);
CFL=0.5;
limiter='superbee';
strang=dt*0.5;

%%
% The main loop
if wb
  wstr=sprintf('Computing %d steps.',N); 	
  wbar=waitbar(0,wstr);
end;
for k=1:N,
  U=central(xflux,U,strang,dx,xdir,boundary_cond,CFL,limiter);
  if wb,
    waitbar((k-0.5)/N,wbar);
  end;
  if (k==1),
    strang=dt;
  end;
  U=central(yflux,U,strang,dy,ydir,boundary_cond,CFL,limiter);
  if k==N,
    strang=0.5*dt;
    U=central(xflux,U,strang,dx,xdir,boundary_cond,CFL,limiter);
  end;
  if wb,
    waitbar(k/N,wbar);
  end;
end;

%%
% 
if wb,
  close(wbar)
end;
rmpath ../../AppendixA;