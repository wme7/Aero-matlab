function u=transheat(u0,a,mu,x,T,nstep)
%    Solves the equation u_t + au_x  = mu u_{xx} by operator
%    splitting.  u0 and a must the size of x.
%
%    Output is a matrix of size lentgth(x)\times (T/dt)+1. 
%    length(x) must be even!
%
%

% Setting up a path to functions from other examples
h=path;
path(h,'../Example2_2');
path(path,'../Example2_3');

% Initial setup
Nt=nstep;
dt=T/Nt;
dx=x(2)-x(1);
m=length(x);
u=zeros(m,Nt+1);

% Solving by splitting
u(:,1)=u0;
for i=1:Nt,
	u1=transport(a,u(:,i),dt,dx,1,'periodic'); % Transport
	u(:,i+1)=heat(u1,mu*dt,1);                 % Diffusion
end;

path(h);