%% Example:
% Compute the tridiagonal matrix used to solve implicitly a laplace (heat)
% equation.

m = 5; % domain size.
h = 0.1; % delta x value.

I = eye(m);
e = ones(1,m);
x = spdiags([e' -4*e' e'],[-1 0 1],m,m);
%x = zeros(m);
%x = x+spdiags([e' -4*e' e'],[-1 0 1],m,m);
s = spdiags([e',e'],[-1 1],m,m);
s = zeros(m);
s = s+spdiags([e',e'],[-1 1],m,m);
%kron(I,x);
%kron(s,I);
A = (kron(I,x)+kron(s,I))/h^2;
%size(A);