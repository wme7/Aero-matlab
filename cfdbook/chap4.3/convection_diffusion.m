% Numerical solution of 1-D convection-diffusion equation with Dirichlet
% boundary conditions:
%	du/dx - (1/Pe)d^2u/dx^2 = 0,  u(0) = a,   u(1) = b
% Central scheme for convection, local grid refinement

% Copyright 2001 P. Wesseling
% This program and its subprograms may be freely used, modified and distributed
% under the GNU General Public License: http://www.gnu.org/copyleft/gpl.html

% Theory is given in Section 4.3 of
% 	P. Wesseling: Principles of Computational Fluid Dynamics
% 	Springer, Heidelberg, 2000 ISBN 3-5453-0. XII, 642 pp.
% See http://ta.twi.tudelft.nl/nw/users/wesseling/cfdbook.html

% This program generates Figures 4.2 and 4.5 in the book

% Functions called: exact_solution

		% .......................Input..............................
a = 0.0;	% Left Dirichlet value
b = 1.0;	% Right Dirichlet value
pe = 40;	% Peclet number
m1 = 8;		% Number of cells outside refinement zone
m2 = 4;		% Number of cells inside refinement zone
refinement = 0;	% Enter 1 for for local grid refinement or something else for
		% uniform grid ..........End of input......................

if refinement == 1
  del = 4/pe;			% del = size of refinement zone
else
  del = m2/(m1+m2);
end

		% Definition of cell numbering.............................

%                              |<---- del----->|
%     x=0       m1 cells       |    m2 cells  x=1
% grid |---o---|---o---|---o---|-o-|-o-|-o-|-o-|
%      1   1   2   2   3   3   4 4 5 5 6 6 7 7 8

x = zeros(m1+m2+1,1);		% x contains coordinates of cell boundaries
for j = 1:(m1+1)
  x(j) = (j-1)*(1-del)/m1;
end
for j = (m1+2):(m1+m2+1),
  x(j) = x(m1+1) + (j-1-m1)*del/m2;
end
n = length(x)-1; dx = zeros(n,1);	% dx contains cell sizes
y = zeros(n,1);				% y contains coordinates of cell centers 
for j = 1:n
  y(j) = 0.5*(x(j)+x(j+1)); dx(j) = x(j+1)-x(j);
end

%	Numerical scheme gives system Au = f
% 	Diagonals of tridiagonal matrix A are stored in B:
% 	A(j,j-1) = B(j-1,1)  A(j,j) = B(j,2)  A(j,j+1) = B(j+1,3)
B = zeros(n,3); f = zeros(n,1);
B(1,2) = 0.5 + (2/pe)*(1/dx(1) + 1/(dx(1)+dx(2)));
B(2,3) = 0.5 - (2/pe)/(dx(1)+dx(2));
f(1)   = a*(1+2/(pe*dx(1)));
for j = 2:(n-1)
  B(j-1,1) = -0.5 - (2/pe)/(dx(j)+dx(j-1));
  B(j+1,3) = 0.5 - (2/pe)/(dx(j)+dx(j+1));
  B(j,2)   = - B(j-1,1) - B(j+1,3);
end
j = n; B(j-1,1) = -0.5 - (2/pe)/(dx(j)+dx(j-1));
       B(j,2)   = -0.5 + (2/pe)*( 1/(dx(j)+dx(j-1)) +1/dx(j));
       f(j)     = -b +(2*b)/(dx(j)*pe);
       
numerical_solution = spdiags(B, [-1 0 1], n,n)\f;

% Computation of error norms
error = numerical_solution - exact_solution(y,pe,a,b); 
n = length(y); norm1 = norm(error,1)/n
norm2 = norm(error,2)/sqrt(n)
norminf = norm(error,inf)

% Make plot
clf, z = 0:0.01:1; plot(z, exact_solution(z,pe,a,b)), hold on
title(['Pe=',num2str(pe),',  ',num2str(n),' cells,  ',num2str(m2),...
' cells in boundary layer'],'fontsize',18)
plot(y,numerical_solution,'*','markersize',13)


