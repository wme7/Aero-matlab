function u=EOBurg(u,t,dx,dir,bndcond)
%%
% Solves Burgers' equation u_t + 0.5*u^2_x = 0 by the Engquist-Osher method.
%

%%
% Initial setup
if nargin<5,
  bndcond='Neumann';
end;
if nargin<4,
  dir=1;
end;
if (dir>2),
  error(strcat('You can easily extend this to',num2str(dir)','D yourself ...'));
end;
% Some things ensuring that you solve along the columns of u. Must
% transpose if u is a row vector or if dir = 2.
transpose=0;
S=size(u);
N=S(1);
if ((N==1)||(dir==2)),
  u=u';
  transpose=1;
end;
maxspeed=max(max(abs(u)));
CFL=0.95;
dt=CFL*dx/maxspeed; 
nstep=ceil(t/dt);
dt=t/nstep;
lambda=dt/dx;
%%
% Solving by EO metod
for j=1:nstep,
  uext=set_extend('iden',u,bndcond);
  u=u-lambda*diff(eoflux(uext));
end;
if (transpose)
  u=u';
end;

%% eoflux
function f=eoflux(u)
  S=size(u);
  N=S(1);
  MM=max(u(1:N-1,:),0.0); mm=min(u(2:N,:),0.0);
  f=0.5*(MM.^2+mm.^2);
  
  
  
  
