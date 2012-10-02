function f=linearfunc(u,deriv)
% f(u)=u to be used with the finite volume code.
if nargin<2,
	deriv=0;
end;
if deriv==0,
	f=u;
else
	f=ones(size(u));
end;