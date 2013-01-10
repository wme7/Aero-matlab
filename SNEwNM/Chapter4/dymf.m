function dyu = dymf(u)
% DYMF Matrix-free partial derivative wrt y;
% homogeneous Dirichlet BC.
%
n2=length(u);
n=sqrt(n2);
h=.5*(n+1);
%
% Turn u into a 2D array with the BCs built in.
%
uu=zeros(n+2,n+2);
vv=zeros(n,n);
vv(:)=u;
uu(2:n+1,2:n+1)=vv;
%
% Compute the partial derivative.
%
dyuu=zeros(n,n);
dyuu=uu(2:n+1,3:n+2)-uu(2:n+1,1:n);
%
% Divide by 2*h.
%
dyu=.5*dyuu(:)/h;
