function s=source(u,kappa)
if nargin<2,
	kappa=0.5;
end;
s=kappa*u.*(1-u).*(u-0.5);
