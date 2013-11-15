function u=sourcesplit(flux,source,u0,x,T,nstep)
%%
%    Solves the equation u_t + flux(u)_x  = source(u) by operator
%    splitting.  u0 and a must the size of x.
%
%    Output is a matrix of size |[lentgth(x),((T/dt)+1)]|. 
%

%% Initial setup
h=path;
path(h,'../Example2_5');
Nt=nstep;
dt=T/Nt;
dx=x(2)-x(1);
m=length(x);
u=zeros(m,Nt+1);
u(:,1)=u0;
%% Solving by operator splitting
for i=1:Nt,
	u1=LxF(flux,u(:,i),dt,dx,1,'periodic');    % Transport solving by LxF
	u(:,i+1)=odewrap(source,u1,dt);            % ODE solving by matlab
end;
path(h);
