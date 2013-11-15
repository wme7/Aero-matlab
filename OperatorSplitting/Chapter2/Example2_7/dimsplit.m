function u=dimsplit(fflux,gflux,u0,x,y,T,nstep)
%%
%  Solves the equation u_t + fflux(u)_x  + gflux(u)_y = 0 by operator
%  splitting.  u0 and a must the size of x.
%
%  Output is a matrix of size [length(x),((T/dt)+1)] . 
%

%% Initial setup
h=path;
path(h,'../Example2_5');
Nt=nstep;
dt=T/Nt;
dx=x(2)-x(1); m=length(x);
dy=y(2)-y(1); n=length(y);
u=zeros(m,n,Nt+1);
u(:,:,1)=u0;

%% Dimensional splitting
%  uses the LxF method in each direction.
for i=1:Nt,
	u1=LxF(fflux,u(:,:,i),dt,dx,1,'periodic');      % Conservation in x-dir
	u(:,:,i+1)=LxF(gflux,u1,dt,dy,2,'periodic');    % Conservation in y-dir
end;

path(h);