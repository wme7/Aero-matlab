%% Dimsplit
function U=DimSplit(U,X,Y,xsolver,ysolver,delta,N,T,varargin)
%   Performs dimensional splitting on the 2D scalar conservation law
%   U_t + xflux(U)_x + yflux(U)_y = 0 with initial data U.
%   Uses Strang splitting with timestep dt=T/N. The domain is the square
%   spanned by X and Y. varargin may contain 'wbar' if waitbar is desired,
%   'periodic' followd by [xmin xmax], [ymin ymax].
%
%% Initial setup
addpath ../../FrontTrackRut; % the home of front tracking...
wb=0; 
per=' ';
xbnd=[];
ybnd=[];
% Checking the optional input
nv=length(varargin);
i=1;
while i<=nv,
	if strcmp(varargin(i),'wbar'),
		wb=1;
		i=i+1;
	elseif strcmp(varargin(i),'periodic'),
		per='periodic';
		i=i+1;
		xbnd=varargin{i};
		i=i+1;
		ybnd=varargin{i};
		i=i+1;
    else
      display('Unknown option in DimSplit');
      i=i+1;
	end;
end;
% setup
dt=T/N;
tol=0.001*delta;  % Tolerance for collision times etc.
[n m]=size(U);
strang=dt*0.5;    % The Strang splitting time 
if wb
  wstr=sprintf('Computing %d steps.',N); 	
  wbar=waitbar(0,wstr);
end;
%% The main loop
for k=1:N,
  for j=1:n,
      %% Solving for each column. 
      % x - direction
    IL=find(abs(diff(U(j,:)))>tol);
    if (~isempty(IL)),
      u=U(j,[IL,length(U(j,:))]);
      x=X(IL+1);
      % For each column, fronttracking
      [u,x]=genFrontTrack(u,x,strang,delta,'Riemann_solver',xsolver,per,xbnd);
      U(j,:)=proj(X,u,x);  % projection to piecewise constant on grid.
    end;
  end;
  if wb,
    waitbar((k-0.5)/N,wbar);
  end;
  if (k==1),
    strang=dt; % No longer in the first step (cf. Strang splitting).
  end;
  for i=1:m,
      %% Solving for each row.
      % y -direction
    IL=find(abs(diff(U(:,i)))>tol);
    if (~isempty(IL)),
      u=U([IL;length(U(:,i))],i)';
      x=Y(IL+1);
      % For each row, fronttracking
      [u,x]=genFrontTrack(u,x,strang,delta,'Riemann_solver',ysolver,per,ybnd);
      U(:,i)=proj(Y,u,x)'; % projection to piecewise constant on grid.
    end;
  end;
  if k==N,      % If last step, use the final half step (cf. Strang splitting).
    strang=0.5*dt;
    for j=1:n,
      % Solving for each column. 
      % x - direction
      IL=find(abs(diff(U(j,:)))>tol);
      if (~isempty(IL)),
        u=U(j,[IL,length(U(j,:))]);
        x=X(IL+1);
        % For each column, fronttracking
        [u,x]=genFrontTrack(u,x,strang,delta,'Riemann_solver',xsolver,xbnd);
        U(j,:)=proj(X,u,x); % projection.
      end;
    end;
  end;
  if wb,
    waitbar(k/N,wbar);
  end;
end;
if wb,
  close(wbar)
end;
rmpath ../../FrontTrackRut;