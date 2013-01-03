% Numerical solution of 1-D convection-diffusion equation 
%	d(uv)/dx - (1/Pe)d^2v/dx^2 = q(x),  v(0) = 0,   v(1) = 1
%       u = 1 - b sin(pi.x)
%       Exact solution: v = a{sin(alpha.pi.x) - sin(alpha.pi)} +
%       		\frac{exp[(x-1)Pe]-exp(-Pe)]}{1-exp(-Pe)}
% Hybrid scheme for convection, local grid refinement
% Cell-centered scheme

% Copyright 2001 P. Wesseling
% This program and its subprograms may be freely used, modified and distributed
% under the GNU General Public License: http://www.gnu.org/copyleft/gpl.html

% Theory is given in Sections 4.4 - 4.6 of
% 	P. Wesseling: Principles of Computational Fluid Dynamics
% 	Springer, Heidelberg, 2000 ISBN 3-5453-0. XII, 642 pp.
% See http://ta.twi.tudelft.nl/nw/users/wesseling/cfdbook.html
% This program solves Exercises 4.5.2 and 4.6.1, and generates Figures 4.6
% and 4.7 in the book

% Functions called: exact_solution, hybrid, rhs, velocity

		% .......................Input..............................
a = 0;%0.2; 	% Coefficients in exact
b = 0;%-0.95; 	% solution and
alpha = 1;	% in velocity function
type = 1;	% Enter 1,2 or 3 for upwind, central or hybrid scheme 
pe = 30;	% Peclet number
m1 = 8;		% Number of cells outside refinement zone
m2 = 4;		% Number of cells inside refinement zone
refinement = 0;	% Enter 1 for for local grid refinement or something else for
		%   uniform grid 
m3 = 1;		% Number of defect corrections
		% ....................End of input......................

if refinement == 1
  del = 6/pe;			% del = size of refinement zone
else
  del = m2/(m1+m2);
end

		% Definition of cell numbering.............................

	%                              |<---- del----->|
	%     x=0       m1 cells       |    m2 cells  x=1
	% grid |---o---|---o---|---o---|-o-|-o-|-o-|-o-|
	%      1   1   2   2   3   3   4 4 5 5 6 6 7 7 8

x = zeros(m1+m2+1,1);
for j = 1:(m1+1)
  x(j) = (j-1)*(1-del)/m1;
end
for j = (m1+2):(m1+m2+1)
  x(j) = x(m1+1) + (j-1-m1)*del/m2;
end
n = length(x)-1; dx = zeros(n,1);	% dx contains cell sizes
y = zeros(n,1);				% y contains coordinates of cell centers 
for j = 1:n
  y(j) = 0.5*(x(j)+x(j+1)); dx(j) = x(j+1)-x(j);
end

u = velocity(x,b);			% Velocity at cell faces
p = zeros(n-1,1);			% Mesh Peclet numbers
for j = 1:n-1
  p(j) = u(j+1)*pe*(dx(j)+dx(j+1))/2;
end
s = hybrid(p,type);			% Switch factors for hybrid scheme

	%	Numerical scheme gives system Av = f
	% 	Diagonals of tridiagonal matrix A are stored in B:
	% 	A(j,j-1) = B(j-1,1)  A(j,j) = B(j,2)  A(j,j+1) = B(j+1,3)

B = zeros(n,3); 
B(1,2) = (2/pe)*(1/dx(1) +  1/(dx(1)+dx(2)) ) + (u(2)+s(1)*abs(u(2)))/2;
B(2,3) =  - (2/pe)/(dx(1)+dx(2)) + (u(2)-s(1)*abs(u(2)))/2;
for j = 2:(n-1)
  B(j-1,1) = -0.5*(u(j)+s(j-1)*abs(u(j)))  - (2/pe)/(dx(j)+dx(j-1));
  B(j+1,3) = 0.5*(u(j+1)-s(j)*abs(u(j+1))) - (2/pe)/(dx(j)+dx(j+1));
  B(j,2)   = (2/pe)/(dx(j)+dx(j-1)) + (2/pe)/(dx(j)+dx(j+1)) ...
             +0.5*(u(j+1)+s(j)*abs(u(j+1))) - 0.5*(u(j)-s(j-1)*abs(u(j)));
end
j = n;
B(j-1,1) =  -0.5*(u(j)+s(j-1)*abs(u(j)))  - (2/pe)/(dx(j)+dx(j-1));
B(j,2)   = - 0.5*(u(j)-s(j-1)*abs(u(j)))  + (2/pe)*( 1/(dx(j)+dx(j-1)) +1/dx(j));
f = rhs(y,a,b,alpha,pe);		% Given right-hand side function
for j = 1:n,
  f(j) = f(j)* dx(j);
end
f(n) = f(n) - u(n+1) + 2/(pe*dx(n));	% Boundary modification
% ERROR: f(n) = f(n) - (1-s(n-1))*u(n+1) + 2/(pe*dx(n));	% Boundary modification
numerical_solution = spdiags(B, [-1 0 1], n,n)\f;

	% Computation of error norms
error = numerical_solution - exact_solution(y,a,alpha,pe);
norm1 = norm(error,1)/n
norm2 = norm(error,2)/sqrt(n)
[norminf,location] = max(abs(error))

% Make plot
clf, z=0:0.01:1; plot(z, exact_solution(z,a,alpha,pe)), hold on
if type == 1
 t = ' upwind';
elseif type == 2
 t = ' central';
else
 t = ' hybrid';
end
title(['Pe=',num2str(pe),t,' scheme ',num2str(n),' cells ',num2str(m3),...
' defcorr'],'fontsize',18)
plot(y,numerical_solution,'*','markersize',10)

if m3 > 0			  % Defect correction
  n = length(y); C = zeros(n,3);  % Matrix generation for central scheme (s=0)
  C(1,2) = (2/pe)*(1/dx(1) +  1/(dx(1)+dx(2)) ) + (u(2))/2;
  C(2,3) =  - (2/pe)/(dx(1)+dx(2)) + (u(2))/2;
  for j = 2:(n-2)
    C(j-1,1) = -0.5*(u(j))  - (2/pe)/(dx(j)+dx(j-1));
    C(j+1,3) = 0.5*(u(j+1)) - (2/pe)/(dx(j)+dx(j+1));
    C(j,2)   = (2/pe)/(dx(j)+dx(j-1)) + (2/pe)/(dx(j)+dx(j+1)) ...
              +0.5*(u(j+1)) - 0.5*(u(j));
  end
  j = n;  C(j-1,1) =  -0.5*(u(j))  - (2/pe)/(dx(j)+dx(j-1));
          C(j,2)   = - 0.5*(u(j))  + (2/pe)*( 1/(dx(j)+dx(j-1)) +1/dx(j));
	  
  for i = 1:m3			% m3 defect corrections
    g = f + ....
    (spdiags(B, [-1 0 1], n,n) - spdiags(C, [-1 0 1], n,n))*numerical_solution;
    numerical_solution = spdiags(B, [-1 0 1], n,n)\g;
  end
  
  plot(y,numerical_solution,'o ','markersize',10)

	  % Computation of error norms
  error = numerical_solution - exact_solution(y,a,alpha,pe);
  norm1 = norm(error,1)/n
  norm2 = norm(error,2)/sqrt(n)
  [norminf,location] = max(abs(error))
end
