function fb=bvpsys(u)
% BVPSYS Two point BVP for two unknown functions.
%        FB=BVPSYS(U) is the residual of the boundary value problem
%        formualted as a first-order system for $U = (v, v')$.
%
%        Problem 7.4, page 187 in
%        Computer Methods for Ordinary Differential 
%        Equations and Differential Algebraic Equations
%        by U. M. Ascher and L. R. Petzold, SIAM 1998.
%
% 
global L
n2=length(u);
fb=zeros(n2,1);
n=n2/2; h=L/(n-1);
f1=zeros(n,1); f2=zeros(n,1);
r=0:n-1; r=r'*h;
%
% separate v and v' from their storage in u
%
v=u(1:2:n2-1); vp=u(2:2:n2);
%
% Set the boundary conditions.
%
f1(1)=vp(1);   % v'(0) = 0
f2(n)=v(n);    % v(L) = 0;
%
f1(2:n)= v(2:n)-v(1:n-1)-h*.5*(vp(2:n)+vp(1:n-1));
%
% v'' = (4/t) v' + (t v  - 1) v
%
% The division by zero really doesn't happen. Fix it up.
%
cof=r; cof(1)=1; cof=4./cof; cof(1)=0;
%
rhs=cof.*vp + (r.*v - 1).*v; 
f2(1:n-1)= vp(2:n)-vp(1:n-1) + h*.5*(rhs(2:n)+rhs(1:n-1));
fb(1:2:n2-1)=f1; 
fb(2:2:n2)=f2;
