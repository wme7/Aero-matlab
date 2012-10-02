function f=burgers(u,deriv)
% f(u)=u^2/2 to be used with the finite volume code.
if nargin<2,
	deriv=0;
end;
if deriv==0,
	f=0.5*u.^2;
else
	f=u;
end;