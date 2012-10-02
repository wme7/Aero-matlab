function f=burger(u,deriv)
if nargin<2,
	deriv=0;
end;
if deriv==0,
	f=0.5*u.^2;
else
	f=u;
end;