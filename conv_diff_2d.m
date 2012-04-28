%% Homework 11
% We wish to solve the following PDE using the approximate factorization
% method:
%
% $\frac{\partial \phi }{\partial t}+U(x,y)\frac{\partial \phi }{\partial
% x}+V(x,y)\frac{\partial \phi }{\partial y}=\alpha \left (
% \frac{\partial^2 \alpha }{\partial x^2}-\frac{\partial^2 \phi }{\partial
% y^2} \right)$
% 
% This is 2D convection-diffusion equation

close all
clear all

%% Grid

ly = [-1,1];                lx =[0,10];
p  = 4;                     q  = 8;
dy = (ly(2)-ly(1))/p;       dx = (lx(2)-lx(1))/q;
y  = ly(1):dy:ly(2);        x  = lx(1):dx:lx(2); % Physical grid
t_end = 1;      dt = 0.1;   n = (t_end-0)/dt;
phi_grid = ones(p+1,q+1,n); % Numerical grid


%% Constants
% assume:
alpha = 2;

%% variables
% assume: $u(x,y)=x*y% and $v(x,y)=x/(y+1)$,
U = y' * x; V = (1./(y+0.1))' * x;

%% Formulate Base matrices
% Ax*
e = ones(q-1,1);
c = 2*ones(q-1,1); c(q-1,1)=-1;
Ax_star = spdiags([e c e],[-1 0 1],q-1,q-1);

% Bx*
e = ones(q-1,1);
c = zeros(q-1,1); c(q-1,1)=1;
Bx_star = spdiags([e c -e],[-1 0 1],q-1,q-1);

% True Ax, Bx, 
I = eye((p-1),(p-1));
Ax = kron(I,Ax_star)/dx^2;
Bx = kron(I,Bx_star)/dx;

% Ay*
e = ones(p-1,1);
c = 2*ones(p-1,1); c(p-1,1)=-1;
Ay_star = spdiags([e c e],[-1 0 1],p-1,p-1);

% By*
e = ones(p-1,1);
c = zeros(p-1,1); c(p-1,1)=-1;
By_star = spdiags([e c -e],[-1 0 1],p-1,p-1);

% True Ay, By
I = eye((q-1),(q-1));
Ay = kron(I,Ay_star)/dy^2;
By = kron(I,By_star)/dy;

% I
I = eye((p-1)*(q-1),(p-1)*(q-1));

% Dx
%U(j,i)=1
Dx = alpha*Ax-Bx;

% Dy
%V(j,i)=1
Dy = alpha*Ay-By;

% f=(I-dt/2*Dx)(I-dt/2*Dy)*phi^n
f = (I+dt/2*Dx)*(I+dt/2*Dy);

% z = (I-dt/2*Dx)\f
%z = (I-dt/2*Dx)\f ;

% phi^n+1 = (I-dt/2*Dy)\z
%(I-dt/2*Dy);

%% Compute $f$ for time n step

%f = (I-dt/2*Dx)*(I-dt/2*Dy);


%% Boundary Conditions


