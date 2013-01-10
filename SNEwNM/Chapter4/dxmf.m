function dxu = dxmf(u)
% DXMF Matrix-free partial derivative wrt x;
% homogeneous Dirichlet BC.
%
n2=length(u);
n=sqrt(n2);
h=1/(n+1);
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
dxuu=zeros(n,n);
dxuu=uu(3:n+2,2:n+1)-uu(1:n,2:n+1);
%
% Divide by 2*h and convert back into a 1D array.
%
dxu=.5*dxuu(:)/h;
