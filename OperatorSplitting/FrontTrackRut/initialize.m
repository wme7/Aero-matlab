%% initialize
function [u,s,t_coll,x]=initialize(u0,x0,T,delta,varargin)
%
%	Initializes the fronts to use by fronttracking. The initial values are 
% 	in the vector u0(1:#components,1:n) and the positions of the 
%   discontinuities in x0(1:n-1). delta gives the accuracy of the Riemann solver
%	and T is the time to front track, this is needed to set the collision time.
%
%   varargin may contain 'Riemann_solver' followed by the name of the
%   solver routine with input (ul,ur,delta) and output [u,s] where with s a
%   vector with the speeds of the fronts and u the states between the
%   fronts. Default 'Riemannsol'.
%   varargin can also contain 'periodic' followed by the interval
%   [xmin,xmax]. 
%   varargin can also contain 'start_time', followed by a real.

%% Initial setup
Periodic=0; t_start=0.0;
Rsolver='Riemannsol';
l=length(varargin);
i=1;
while i<=l,
	if strcmp(varargin(i),'periodic')
		Periodic=1;
		i=i+1;
		bndr=varargin(i);
		xleft=bndr{1};
		i=i+1;
		bndr=varargin(i);
		xright=bndr{1};
	elseif strcmp(varargin(i),'Riemann_solver')
		i=i+1;
		Rsolver=varargin{i};
	elseif strcmp(varargin(i),'start_time')
		i=i+1;
		t_start=varargin{i};
	else
		error(' Unknown option in initialize');
	end;
	i=i+1;
end;
[nc n]=size(u0);

%% Solving the Riemann problems
i=1;ip=2;
if n>2,
  while (max(abs(u0(:,i)-u0(:,ip)))<delta*0.001)&&(ip<n),
    i=ip;
    ip=ip+1;
  end;
end;
ni=i;
if (max(abs(u0(:,ni)-u0(:,ni+1)))<=delta*0.001),
  u=u0;
  x=x0;
  t_coll=(T+1)*ones(length(x)+1);
  s=0.0*ones(size(x));
  return;
else
  [ur s]=feval(Rsolver,u0(:,ni),u0(:,ni+1),delta);  % The Riemann solution
                                                    % with states ur and 
                                                    % speeds s.
  while isempty(s),
	  ni=ni+1;
	  if ni==n,
		  u=u0;
		  x=x0;
		  t_coll=(T+1)*ones(length(x)+1);
		  s=0.0*ones(size(x));
		  return;
	  end;
	  [ur s]=feval(Rsolver,u0(:,ni),u0(:,ni+1),delta);
  end;
  u=[u0(:,ni) ur u0(:,ni+1)];
  [nc nlast]=size(u);
  x=x0(ni)*ones(size(s));
  t_coll=(T+1)*ones([1,nlast]);
  for i=ni+1:n-1
    [ur sri]=feval(Rsolver,u(:,nlast),u0(:,i+1),delta);
    nx=length(x);
    if ~isempty(sri),
      u=[u ur u0(:,i+1)];
      s=[s sri];
      % Finds the collision time.
      t_coll(nlast)=collision_time([s(nlast-1) sri(1)],[t_start,t_start],[x(nx),x0(i)]);
      ns=length(sri);
      t_coll=[t_coll (T+1)*ones([1,ns])];
      x=[x x0(i)*ones([1,ns])];
    end;
    [nc nlast]=size(u);
  end;
  if Periodic, % Special considerations if the boundary conditions are periodic.
	   if x(1)<=xleft,
		  error('  Not proper left boundary for periodic');
	  end;
	  if x(nlast-1)>=xright,
		  error(' Not proper right boundary for periodic');
	  end;
	  [ur sri]=feval(Rsolver,u(:,nlast),u(:,1),delta);
	  if ~isempty(sri),
		  ns=length(sri);
		  i=1;
		  UR=[u(:,nlast) ur u(:,1)];
		  if sri(ns)<=0.0,
			  i=ns+1;
		  else
			  while (sri(i)<=0.0),
				  i=i+1;
			  end;
		  end;
		  UN=[UR(:,2:i) UR(:,i)];
		  U1=[UR(:,i) UR(:,i:ns)];
		  sN=[sri(1:i-1) 0.0];
		  s1=[0.0 sri(i:ns)];
		  tcoll_left=collision_time([s1(length(s1)) s(1)],...
			  [t_start t_start],[xleft x(1)],T);
		  tcoll_right=collision_time([s(nlast-1) sN(1)],...
			  [t_start t_start], [x(nlast-1) xright],T);
		  n1=ones(size(s1));
		  nN=ones(size(sN));
		  x=[xleft*n1 x xright*nN];
		  s=[s1 s sN];
		  
		  t_coll=[(T+1)*n1 tcoll_left ...
			  t_coll(2:nlast-1)...
			  tcoll_right (T+1)*nN];
		  u=[U1 u UN];
	  else
		  tcoll_left=collision_time([0 s(1)],[t_start t_start],[xleft x(1)],T);
		  tcoll_right=collision_time([s(nlast-1) 0.0],...
			  [t_start t_start],[x(nlast-1) xright],T);
		  s=[0 s 0];
		  x=[xleft x xright];
		  t_coll=[T+1 tcoll_left t_coll(2:nlast-1) tcoll_right T+1];
		  u=[u(:,1) u u(:,nlast)];
	  end;
  end;
end;

