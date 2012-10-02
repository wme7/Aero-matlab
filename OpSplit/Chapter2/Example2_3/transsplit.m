function u=transsplit(u0,a,b,x,y,T,nstep)
%%
% Computes an approximation to the solution of 
%
% $$ u_t+  a u_x + b u_y=0 $$
%
% with initial data u0 on the rectangle spanned by the vectors x and y. u0 must be
% of the size |[length(x),length(y)]|.
%
% Output is a 3d matrix of size |[lentgth(x),length(y),((T/dt)+1)]|.
%
%% Initial setup
dt=T/nstep; Nt=ceil(T/dt); dt=T/Nt;
dx=x(2)-x(1); dy=y(2)-y(1);
u=zeros(length(x),length(y),Nt+1);

%% Solving by Strang splitting 
u(:,:,1) = u0;
u1 = transport(b,u(:,:,1),0.5*dt,dx,1);   % Transport in the x-direction
for i=2:Nt,
   u(:,:,i) = transport(a,u1,      dt,dy,2);   % Transport in the y-direction
   u1       = transport(b,u(:,:,i),dt,dx,1);   % Transport in the x-direction
end
u(:,:,Nt+1) = transport(a,u1,            dt,dy,2);   % Transport in the y-direction
u(:,:,Nt+1) = transport(b,u(:,:,Nt+1),0.5*dt,dx,1);   % Transport in the x-direction

