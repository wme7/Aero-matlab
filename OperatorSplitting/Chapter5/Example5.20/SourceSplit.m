%% Sourcesplit
function U=SourceSplit(u0,X,solver,sourcefnc,delta,N,T,varargin)

%   Performs operator splitting on the scalar conservation law
%   U_t + flux(U)_x  = sourcefnc(u) with initial data u0. The two 
%   equations are U_t + flux(U)_x = 0, solved with front tracking, and 
%   U_t = source(U), solved with Euler's method. The Riemann solver is
%   specified by 'solver'. 
%   Uses Strang splitting with timestep dt=T/N. The domain is the line
%   spanned by X. varargin may contain 'wbar' if waitbar is desired,
%   'periodic' followd by [xmin xmax], 'central' if a central method is to
%   be used instead of front tracking, and 'kappa' followed by a parameter
%   to the sourcefnc.
%    

%% Initial setup
addpath ../../FrontTrackRut;
addpath ../../AppendixA
wb=0; 
per=' ';
xbnd=[];
Xf=X;
nv=length(varargin);
heun=0;
i=1;
fronttrack=1;
strangsplit=0;
NT=0;
while i<=nv,     % reading varargin
	if strcmp(varargin(i),'wbar'),
		wb=1;
		i=i+1;
	elseif strcmp(varargin(i),'heun'),
		heun=1;
		i=i+1;
	elseif strcmp(varargin(i),'periodic'),
		per='periodic';
		i=i+1;
		xbnd=varargin{i};
		i=i+1;
		nx=length(X);
		Xf=X(2:nx-1);
	elseif strcmp(varargin(i),'kappa'),
		i=i+1;
		Kappa=varargin{i};
		i=i+1;
	elseif strcmp(varargin(i),'fronttrack'),
		i=i+1;
	elseif strcmp(varargin(i),'central'),
		fronttrack=0;
		NT=1;
		i=i+1;
	elseif strcmp(varargin(i),'strang'),
		strangsplit=1;
		i=i+1;
	else
      display('Unknown option in SourceSplit');
	  i=i+1;
	end;
end;

dt=T/N;
n=length(u0);
U=zeros(N+1,n);
U(1,:)=u0;

%% Solving
if strangsplit,
	strang=dt*0.5;
else
	strang=dt;
end;
if wb
  wstr=sprintf('Computing %d steps.',N); 	
  wbar=waitbar(0,wstr);
end;
for k=1:N,
	u=U(k,:);
	x=Xf;
	%Front-tracking
	if fronttrack,
		[u,x]=genFrontTrack(u,x,strang,delta,'Riemann_solver',solver,per,xbnd);
        % Front tracking
		w=proj(X,u,x);  % projection.
	elseif NT,
		dx=x(2)-x(1);
		w=central('burgers',u,strang,dx,1,per,0.5,'superbee');
        % Using a central scheme
	end;
	if wb,
		waitbar((k-0.5)/N,wbar);
	end;
	if (k==1),
		strang=dt;
	end;
	fw=feval(sourcefnc,w,Kappa);
	U(k+1,:)=w+dt*fw;%  Euler step for sourceterm
	if heun,
		U(k+1,:)=w+0.5*dt*(feval(sourcefnc,U(k+1,:),Kappa)+fw);  % Heun step
	end;
	if strangsplit,
		if k==N,
			strang=0.5*dt;
			u=U(k+1,:);
			x=Xf;
			% Front-tracking
			if fronttrack,
				[u,x]=genFrontTrack(u,x,strang,delta,'Riemann_solver',solver,per,xbnd);
				U(k+1,:)=proj(X,u,x); % projection.
			elseif NT,
				U(k+1,:)=central('burgers',u,strang,dx,1,per,0.5,'mm');
			end;
		end;
	end;
	if wb,
		waitbar(k/N,wbar);
	end;
end;

%% Closing ...
if wb,
  close(wbar)
end;
rmpath ../../FrontTrackRut;
rmpath ../../AppendixA