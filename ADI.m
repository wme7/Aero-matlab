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
dy = 0.2;                   dx = 0.5;
y  = ly(1):dy:ly(2);        x  = lx(1):dx:lx(2); % Physical grid
q  = length(y);             p  = length(x); 
phi = zeros(q,p); % Numerical grid

t_end = 1;      dt = 0.1;   n = (t_end-0)/dt;

%% Boundary Conditions

phi(:,1) = 1;
phi(1,:) = 1;


%% Constants & coheficients
% assume: $\alpha = const$.

alpha = 1;

% 1. assume: $u(x,y)=x*y% and $v(x,y)=x/(y+0.1)$,
% 2. assume: $u(x,y)=1$ and $v(x,y)=1$.

u = zeros(q,p); 
for j = 1:q
    for i = 1:p
        u(j,i) = x(i)*y(j);
        % u(j,i) = 1 
    end
end

v = zeros(q,p);
for j = 1:q
    for i = 1:p
        v(j,i) = (1/(y(j)+0.1))*x(i);
        % v(j,i) = 1 
    end
end

%% Formulate Base matrices
% Dy
% Constants
a = dt/2*(alpha/dy^2 - v/(2*dy));
b = (1 - dt/dy^2);
c = dt/2*(alpha/dy^2 + v/(2*dy));

Qy = zeros(q-2,q-2,p);
for i = 2:p-1
    e = zeros(q-2,3);
    e(:,3) =  a(1:q-2,i);
    e(:,2) =  b*ones(q-2,1); % But the last No. in de Diag must be modif.:
    e(q-2,2) =  a(q-1,i)-b(q-1,i);
    e(:,1) =  c(3:q,i);
    Qy(:,:,i) = spdiags(e,-1:1,q-2,q-2);
end


% Dx
% Constants
% f = dt/2*(alpha/dx^2 - u/(2*dx));
% g = (1 - dt/dx^2);
% h = dt/2*(alpha/dx^2 - u/(2*dx));
% 
% e1 = a(1,3:p-1);
% e2 = b*ones(1,p-2);
% e3 = c(1,2:p-2);
% Qx = spdiags([e3 e2 e1],[-1 0 1],p-2,p-2);



%% Compute $f$ for time n step



%% Boundary Conditions


