function f=fflux(u,i)
u2=u.^2;
n=u2+(1-u).^2;
if nargin<2 || (i==0)
	% give the flux 
	f= u2./n;
elseif (i==1)
    % the derivative 
	f= 2*u.*(1-u)./(n.^2);
else
	error(' Unknown option in fflux ');
end;