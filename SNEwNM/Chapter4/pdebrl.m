% PDEBRL Solve the left preconditionied PDE with brsola.
% This program creates the left preconditionied PDE example
% from Chapter 4
%
%
% The example has the form
%
% L(u) = f, 0 < x,y < 1
% with homogeneous Dirichlet BC
%
% In terms of the program below, the equation is
% Lu = - (u_xx + u_yy) + 20*u*(u_x + u_y)  = f
%
% Within the nonlinear map, pdefun.m, we apply fish2d.m as
% a preconditioner.
%
% The action of the laplacian, partial wrt x, and partial wrt y
% are computed in the routines lapmf, dxmf, and dymf. Here mf
% stands for MATRIX FREE.
%
global rhsf prhsf;
m=31; n=m*m;
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
% Form the solution: sol(x,y)= 10 x y (1-x) (1-y) exp(x^4.5)
%
zsoly=x.*(1-x);
zsolx=zsoly.*exp(x.^4.5);
solt=10*zsolx*zsoly';
sol=solt(:);
clear rhs; clear axt; clear at; clear ayt; clear solt;

% fix the left side so that the solution is known
%
% unpreconditioned and left preconditioned
%
rhsf = lapmf(sol) + 20*sol.*(dxmf(sol)+dymf(sol));
%
% left preconditioned
%
prhsf = sol + 20*fish2d(sol.*(dxmf(sol)+dymf(sol)));
%
% Set the initial iterate to zero. Call the nsoli solver with GMRES
% and the compare to Broyden's method
%
u0=zeros(m*m,1);
%
tol=[1.d-1,1.d-1]*tolh;
%
% GMRES (default parameters)
%
x=zeros(n,1); 
[sol, it_histg, ierr] = nsoli(x,'pdeleft',tol);
%
% BiCGSTAB
%
parms = [40,40,.9,3];
x=zeros(n,1); 
parms=[40,40,-.1];
[sol, it_histg2, ierr] = nsoli(x,'pdeleft',tol,parms);
%
% Broyden's Method
%
x=zeros(n,1); 
[sol, it_histb, ierr] = brsola(x,'pdeleft',tol);
figure(1)
ng=length(it_histg(:,1));
ng2=length(it_histg2(:,1));
nb=length(it_histb(:,1));
semilogy(0:ng-1,it_histg(:,1)/it_histg(1,1),'-',...
0:ng2-1,it_histg2(:,1)/it_histg2(1,1),'-.',...
0:nb-1,it_histb(:,1)/it_histg(1,1),'--');
xlabel('Nonlinear iterations');
ylabel('Relative residual');
legend('GMRES','BICGSTAB','Broyden')
figure(2)
semilogy(it_histg(:,2),it_histg(:,1)/it_histg(1,1),'-',...
it_histg2(:,2),it_histg2(:,1)/it_histg2(1,1),'-.',...
it_histb(:,2),it_histb(:,1)/it_histg(1,1),'--');
ylabel('Relative residual');
xlabel('Function Evaluations');
legend('GMRES','BICGSTAB','Broyden')

