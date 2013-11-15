function u=diffBurg(Afunc,u0,nu,dx,T,Nt,bndcond)
%
%    Solves the equation u_t + 0.5(u^2)_x  = \nu A(u)_{xx} by operator
%    splitting.  u0 and a must the size of x.
%
%    Output is a matrix of size lentgth(x)\times (T/dt)+1. length of x
%    must be even!
%
%%
% initial setup
if nargin<7
	bndcond='neumann';
end;
dt=T/Nt;
m=length(u0);
u=zeros(m,Nt+1);
u(:,1)=u0;

%% 
% Solving by operator splitting
for i=1:Nt,
	u1=EOBurg(u(:,i),dt,dx,1,bndcond);                % Burgers' equation
	u(:,i+1)=nonlindiff(Afunc,nu,u1,dt,dx,1,bndcond); % Nonlinear diffusion
end;

