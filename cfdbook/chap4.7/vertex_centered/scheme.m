% Numerical solution of 2-D singular perturbation problem for stationary
%  convection-diffusion equation:
%  	- dw/dx - D*Laplacian(w) = q
%	x = 1: Dirichlet   x = 0: homogeneous or inhomogeneous Neumann
%	y = 0: Dirichlet   y = 1: homogeneous Neumann
%
% Central or upwind scheme for convection
% Vertex-centered grid; uniform in x, locally refined in y

% Copyright 2001 P. Wesseling
% This program and its subprograms may be freely used, modified and distributed
% under the GNU General Public License: http://www.gnu.org/copyleft/gpl.html

% Theory is given in Section 4.7 of
% 	P. Wesseling: Principles of Computational Fluid Dynamics
% 	Springer, Heidelberg, 2000 ISBN 3-5453-0. XII, 642 pp.
% See http://ta.twi.tudelft.nl/nw/users/wesseling/cfdbook.html

% This program generates tables 4.5 and 4.6 in the book

% Functions called: exact_solution, fL, q

		% ...................Input...............................
D = 0.001;	% Diffusion coefficient ( = 1/[Peclet number] )
nx = 8;		% Horizontal number of cells
nyc = 4;	% Vertical number of cells outside refinement zone
nyf = 32;	% Vertical number of cells inside refinement zone
central = 0;	% Enter 1 for central scheme or something else for upwind scheme
hom = 0;	% Neumann condition at outflow boundary: 1 for homogeneous,
		% or something else for inhomogeneous (exact) condition
th = 8;		% Governs thickness of refinement zone = del = th*sqrt(D)
		% .....................End of input......................

del = th*sqrt(D);	% Thickness of refinement zone

% Implementation of upwind  scheme by artificial diffusion in x-direction
dx = 1/(nx-1);
if central == 1
  Da = D;		% Artificially increased diffusion coefficient
else
  Da = D + dx/2;
end

% 		Definition of volume boundaries and grid points
 
%     x=0                                     x=1
% grid o---|---o---|---o---|---o---|---o---|---o
%      1       2                               J

%      |<---- del----->|                      
%     y=0    K1 vert.  |      K2 vertices     y=1
% grid o-|-o-|-o-|-o-|-o---|---o---|---o---|---o
%      1   2   3       K1                      K

J = nx; K1 = nyf; K2 = nyc;
x = [0:J-1]'*dx;	% Horizontal coordinates of grid points (cell vertices) 

K = K1 + K2;
y = zeros(K,1);		% Vertical coordinates of grid points (cell vertices)
ddy = del/(K1-1);
for k = 1:K1,  y(k) = (k-1)*ddy; end
ddy = (1-del)/K2;
for k = K1+1:K,  y(k) = del + (k-K1)*ddy; end
dy = zeros(K,1);	% Vertical distances between vertices
dy(1) = del/(K1-1);
for k = 2:K,  dy(k) = y(k) - y(k-1); end

rhs = zeros(J*K,1);	% Right-hand side of differential equation
for k = 2:K-1		% Lexicographical ordering   of unknowns:
  for j = 1:J-1 	% 	j: x-direction  k: y-direction
    rhs(j + (k-1)*J) = q(x(j),y(k),D)*dx*(dy(k)+dy(k+1))*0.5;
  end
end
k = K;  for j = 1:J-1, rhs(j + (k-1)*J) = q(x(j),y(k),D)*dx*0.5*dy(k); end

for k = 2:K		% Modification of rhs due to halved volumes at x = 0
  m = 1 + (k-1)*J;
  rhs(m) = 0.5*rhs(m);
end
	% Modification of rhs due to boundary conditions
	% x = 0: outflow: homogeneous or inhomogeneous Neumann, depending on 
	% 	 parameter hom
for k = 2:K-1
  m = 1 + (k-1)*J; rhs(m) = rhs(m) + 0.5*(dy(k)+dy(k+1))*Da*fL(y(k),D,hom);
end
k = K;  m = 1 + (k-1)*J; rhs(m) = rhs(m) + 0.5*dy(k)*Da*fL(y(k),D,hom);

	% x = 1: Dirichlet
for k = 1:K, rhs(k*J) = exact_solution(1,y(k),D); end

	% y = 1: homogeneous Neumann,  y = 0: Dirichlet
for j = 1:J,  rhs(j) = exact_solution(x(j),0,D); end

	% Matrix construction by flux evaluation
n = J*K; z = zeros(n,1); d = [-J  -1  0  1  J];
A = spdiags([ z  z  z  z  z], d',n,n);	   % Declaration of sparsity pattern

	% x = 1: Dirichlet
for k = 1:K,  m = k*J;  A(m,m) = 1; end
	% y = 0: Dirichlet
for j = 1:J,  A(j,j) = 1; end

for k = 2:K-1			% Loop over vertical faces
  m1 = (k-1)*J;  m = 1 + m1;  A(m,m) = A(m,m) + 0.5*(dy(k) + dy(k+1));
  a = - (dy(k) + dy(k+1))*(0.25 + 0.5*Da/dx);
  b = - (dy(k) + dy(k+1))*(0.25 - 0.5*Da/dx);
  for j = 2:J-1
    m = j + m1; A(m,m)     =  A(m,m)     - a; A(m,m-1) =  A(m,m-1) - b;
                A(m-1,m-1) =  A(m-1,m-1) + b; A(m-1,m) =  A(m-1,m) + a;
  end
  m = J+m1;  A(m-1,m-1) =  A(m-1,m-1) + b;
	     A(m-1,m)   =  A(m-1,m)   + a;
end
k = K;  m1 = (k-1)*J;  m = 1 + m1;  A(m,m) = A(m,m) + 0.5*dy(k);
  a = - dy(k)*(0.25 + 0.5*Da/dx);
  b = - dy(k)*(0.25 - 0.5*Da/dx);
  for j = 2:J-1			% Loop over horizontal faces
    m = j + m1;
    A(m,m)     =  A(m,m) - a;        A(m,m-1) =  A(m,m-1) - b;
    A(m-1,m-1) =  A(m-1,m-1) + b;    A(m-1,m) =  A(m-1,m) + a;
  end
  m = J+m1;  A(m-1,m-1) =  A(m-1,m-1) + b;  A(m-1,m) =  A(m-1,m) + a;

j = 1; k = 2;			% Loop over horizontal faces;
    m = j + (k-1)*J;                a = - 0.5*dx*D/dy(k-1);
    A(m,m) = A(m,m) - a;            A(m,m-J) = A(m,m-J) + a;
  for k = 3:K
    m = j + (k-1)*J;                a = - 0.5*dx*D/dy(k-1);
    A(m,m) = A(m,m) - a;            A(m,m-J) = A(m,m-J) + a;
    A(m-J,m-J) = A(m-J,m-J) - a;    A(m-J,m) = A(m-J,m) + a;
  end
for j = 2:J-1
k = 2;
    m = j + (k-1)*J;    	    a = - dx*D/dy(k-1);
    A(m,m) = A(m,m) - a;            A(m,m-J) = A(m,m-J) + a;
  for k = 3:K
    m = j + (k-1)*J;                a = - dx*D/dy(k-1);
    A(m,m) = A(m,m) - a;            A(m,m-J) = A(m,m-J) + a;
    A(m-J,m-J) = A(m-J,m-J) - a;    A(m-J,m) = A(m-J,m) + a;
  end
end

numerical_solution = A\rhs;

exsol = zeros (K*J,1);			% Exact solution
for k = 1:K
  for j = 1:J,    exsol(j + (k-1)*J) = exact_solution(x(j),y(k),D);  end
end

% Text output
error = numerical_solution - exsol;  n = K*J;
norm1 = norm(error,1)/n
norm2 = norm(error,2)/sqrt(n)
[norminf,m] = max(abs(error))

	% Solution profile at x = x(1)
hh = zeros(K,1);
for k = 1:K,  hh(k) = numerical_solution(1 + (k-1)*J); end
figure (1), clf, subplot(2,2,1), hold on, plot(hh,y,'o')
yy = 0:0.02:1; exact = zeros (length(yy));
for kk = 1:length(yy),  exact(kk) = exact_solution(x(1), yy(kk), D); end
plot(exact,yy)

	% Solution profile at y = y(kkk)
vv = zeros(J,1); kkk = 1;
for j = 1:J,  vv(j) = numerical_solution(j + (kkk-1)*J); end
subplot(2,2,4), plot(x,vv,'o')
for kk = 1:length(yy),  exact(kk) = exact_solution( yy(kk), y(kkk), D); end
hold on, plot(yy,exact)

% Solution profile at x = x(J)
hh = zeros(K,1); for k = 1:K,  hh(k) = numerical_solution(J + (k-1)*J); end
subplot(2,2,2), plot(hh,y,'o'), hold on
for kk = 1:length(yy),  exact(kk) = exact_solution(x(J), yy(kk), D); end
plot(exact,yy)

	% Contour plot
consol = zeros(K,J);
for k = 1:K
  m1 = (k-1)*J;
  for j = 1:J,    consol(k,j) = numerical_solution(m1 + j);  end
end
subplot(2,2,3), contour(x,y,consol)

numsol = zeros(K,J);		% waterfall plot
for k = 1:K,for j = 1:J, numsol(k,j) = numerical_solution(j + (k-1)*J);end,end
figure(2),clf, waterfall(x,y,numsol), view(200,70), colormap([0 0 0]), grid off






