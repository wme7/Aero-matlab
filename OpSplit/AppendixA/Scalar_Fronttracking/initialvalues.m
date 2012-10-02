function [u,x]=initialvalues(uflux,a,b,init_f,N);
%
%   Approximates the function init_f on the interval [a b] by a piecewise
%   constant function taking values in the set {uflux}.
%   The output is the location of the discontinuities x and u such that 
%   the discontinuity between u(i-1) and u(i) is located at x(i).
%  
if nargin<5,
  N=1/(uflux(2)-uflux(1));
end;
x=linspace(a,b,N);
xf=0.5*(x(2:N)+x(1:N-1));
u1=feval(init_f,xf);
A=ones(size(u1))'*uflux;
B=ones(size(uflux))'*u1;
[m ind]=min(abs(A'-B));
h=uflux(ind);
x=x(2:N-1);
d=diff(h);
ind=d~=0;
[i s]=find(ind);
u=[h(s) h(N-1)];
x=x(s);

