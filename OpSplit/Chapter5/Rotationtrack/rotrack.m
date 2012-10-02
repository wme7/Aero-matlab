function u=rotrack(u0,x,y,nu,T) 
%
%  Solves the linear equation u_t + y u_x - x u_y = 0 by "front tracking" +
%  dimensional splitting, using a time step dictated by dt/dx=nu for each
%  substep. 
%  Neumann boundary condition are implicit. 
%   
dx=x(2)-x(1);
dt=nu*dx;
nstep=ceil(T/dt);
dt=T/nstep;
nu=dt/dx;
u=u0;
for i=1:nstep,
	u=transp(u,nu*y,1);
	u=transp(u,-nu*x,2);
end;

