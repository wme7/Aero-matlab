function f=burger(u,der)
% if second arg. absent u^2/2, else u.
if nargin<2
	f=0.5*u.^2;
else
	f=u;
end;