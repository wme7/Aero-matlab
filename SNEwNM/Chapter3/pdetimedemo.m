% PDETIMEDEMO
% This program creates the left preconditionied time-dependent PDE example
% from Chapter 3
%
%
% The example has the form
%
% u_t = - L(u) + f, 0 < x,y < 1
% with homogeneous Dirichlet BC
%
% In terms of the program below, the equation is
% Lu = - (u_xx + u_yy) + 20*u*(u_x + u_y) 
%
% Within the nonlinear map, pdetime.m, we apply fish2d.m as
% a preconditioner. 
%
% The action of the Laplacian, partial wrt x, and partial wrt y
% are computed in the routines lapmf, dxmf, and dymf. Here mf
% stands for MATRIX FREE.
%
global rhsf uold dt;
m=63; n=m*m;
h=1/(m+1);
h2=(m+1)*(m+1);
h1=(m+1);
tolh=1/h2;
mres=3;
%
% set up the equation
%
x=1:m;
e=x';
x=x'*h;
%
% Form the steady-state solution: sol(x,y)= 10 x y (1-x) (1-y) exp(x^4.5)
%
zsoly=x.*(1-x);
zsolx=zsoly.*exp(x.^4.5);
solt=10*zsolx*zsoly';
sol=solt(:);
clear rhs; clear axt; clear at; clear ayt; clear solt;
%
% Fix the right side so that the steady-state solution is known
%
rhsf = lapmf(sol) + 20*sol.*(dxmf(sol)+dymf(sol));
%
% Set the initial data to zero. 
%
u0=zeros(m*m,1);
uold=u0;
uinit=uold;
dt=.1;
nt=10;
%
%
% Use nsoli.m with GMRES as the solver. 
% Converged result from previous time step is initial iterate.
%
for it=1:nt+1
    norm(uold-sol)
    tol=[1.d-1,1.d-1]*tolh;
    x=zeros(n,1); parms = [40,40,-.1,1];
    [unew, it_histg, ierr] = nsoli(uinit,'pdetime',tol,parms);
    uinit=uold;
    uold=unew;
end
norm(uold-sol)

