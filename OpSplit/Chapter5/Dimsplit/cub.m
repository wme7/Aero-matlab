function f=cub(u,der)
% if second arg. absent then u^3/3 else u^2.
if nargin<2,
	f=0.333333*u.^3;
else
	f=u.^2;
end;