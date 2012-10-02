function u=nonlindiff(Afunc,epsilon,u,t,dx,dir,bndcond)
%%
% Approximating solutions of u_t = \nu A(u)_{xx} by a Cranck-Nicholson
% scheme.  The initial value is u, Afunc is the name of A, must have
% syntax, y=Afunc(u,i), where output is the flux if "i" is zero or
% absent, and the derivative if "i" is one. 
%

%% 
% Initial setup
if (dir>2),
  error(strcat('You can easily extend this to ',num2str(dir),'D yourself ...'));
end;

if nargin<7,
  bndcond='Neumann';
end;
if nargin<6,
  dir=1;
end;
% We are solving along the columns of u
transpose=0;
S=size(u);
N=S(1);
if (N==1||dir==2),
  u=u';
  transpose=1;
end;

%%
% ------------------------------------------------------------------
%
% Approximating the solution  by the C-N method, here "theta" is a
% parameter in [0 1] giving the "implicitness", theta=0 gives an explicit
% method and theta=1 gives an implicit method. 
%
theta=0.5;

speed=feval(Afunc,u,1);
if min(speed)<0,
  error(' Afunc must be nondecreasing for well posedness.');
end;
maxspeed=max(max(speed));
CFL=0.23;                        % Since we are using simple iteration to solve
                                 % the nonlinear part, we must have
                                 % CFL<1/4.
dt=CFL*dx*dx/(theta*maxspeed*epsilon); 
nstep=ceil(t/dt);
dt=t/nstep;
mu=epsilon*dt/(dx.^2);
for j=1:nstep,
  aext=set_extend(Afunc,u,bndcond);
  rhs=u+mu*(1-theta)*diff(diff(aext));
  u=nonlinsolve(Afunc,theta*mu,rhs,u,bndcond);
end;
if (transpose)
  u=u';
end;
	
%-----------------------------------------------------------
%-----------------------------------------------------------

%% nonlinsolve 
function w=nonlinsolve(Afunc,par,rhs,w0,bndcond,tol)
%
% Solving w-par*diff(diff(Afunc(w)))=rhs by simple iteration.
% Terminates if the difference between sucessive iterations is less than
% tol  
  
if nargin<6,
  tol=1e-5;
end;
maxiter=120;
w1=rhs+par*diff(diff(set_extend(Afunc,w0,bndcond)));
error=sum(sum(abs(w1-w0)));
iter=1;
while (error>tol)&&(iter<maxiter),
  w0=w1;
  w1=rhs+par*diff(diff(set_extend(Afunc,w0,bndcond)));
  error=sum(sum(abs(w1-w0)));
  iter=iter+1;
end;
if (iter==maxiter)
  warning('Possibly not convergence in nonlinsolve %10.7f\n',error);
end;
w=w1;
