function u=consheat(flux,u0,x,T,nstep, epsilon)
%%
%    Solves the equation u_t + flux(u)_x  = \eps u_{xx} by operator
%    splitting.  u0 and flux must equal the size of x.
%
%    Output is a matrix of size lentgth(x)\times (T/dt)+1. length of x
%    must be even!
%

%%
% Initial setup
if nargin<6,
   epsilon=0.1; % A parameter for the size of the diffusion
end

% Adding path to previous examples
h=path;
path(h,'../Example2_2');
Nt=nstep;
dt=T/Nt;
dx=x(2)-x(1);
m=length(x);
u=zeros(m,Nt+1);
u(:,1)=u0;

%% 
% Solving by a Lax-Friedrichs method and a spectral diffusion method.
for i=1:Nt,
	u1=LxF(flux,u(:,i),dt,dx,1,'periodic');         % Conservation
	u(:,i+1)=heat(u1,epsilon*dt,1);                 % Small Diffusion
end;

path(h);